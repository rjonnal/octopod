import sys,os
import logging
import numpy as np
import scipy as sp
from scipy.ndimage.filters import generic_filter
from scipy.interpolate import bisplrep,bisplev,interp2d
from matplotlib import pyplot as plt
from scipy.ndimage.morphology import grey_opening
from scipy.ndimage.filters import median_filter
from utils import translation,translation1,autotrim_bscan,find_peaks,shear,Clock,lateral_smooth_3d,polyfit2d,polyval2d,gaussian_convolve
from octopod.Misc import H5
import octopod_config as ocfg
import logging
logging.basicConfig(level=logging.DEBUG)
import time

class BScanMaker:

    def __init__(self,h5,volume_index=0,factor=(1,1),data_block='processed_data',debug=False):
        """Initialize bscan maker. May pass an h5py.File object or .hdf5 filename.
        The file or object must contain a 'processed_data' dataset containing
        at least one volume."""
        self.logger = logging.getLogger(__name__)
        if type(h5)==str:
            self.h5 = H5(h5)
            self.logger.info('Opening file %s.'%self.h5)
        else:
            self.h5 = h5
        self.volume = np.abs(self.h5[data_block][volume_index,:,:,:])
        self.isos = self.h5['projections']['ISOS'][volume_index,:,:]
        self.cost = self.h5['projections']['COST'][volume_index,:,:]
        self.cones = (self.isos+self.cost)/2.0

        # crop self.cones and self.volume
        x,y = self.click_collector(self.cones,'Click upper left and lower right corners for cropping')
        x1 = np.min(x)
        x2 = np.max(x)
        y1 = np.min(y)
        y2 = np.max(y)

        sy,sx = self.cones.shape
        y1 = max(0,y1)
        y2 = min(sy,y2)
        x1 = max(0,x1)
        x2 = min(sx,x2)
        
        self.isos = self.isos[y1:y2,x1:x2]
        self.cost = self.cost[y1:y2,x1:x2]
        self.cones = self.cones[y1:y2,x1:x2]
        self.volume = self.volume[y1:y2,:,x1:x2]


        # select region for B-projection
        x,y = self.click_collector(self.cones,'Click top and bottom edges of region for B-scan')
        x1 = np.min(x)
        x2 = np.max(x)
        y1 = np.min(y)
        y2 = np.max(y)
        sy,sx = self.cones.shape
        y1 = max(0,y1)
        y2 = min(sy,y2)

        subvol = self.volume[y1:y2,:,:]

        subvol_scaled = []
        profs = []
        for iy in range(subvol.shape[0]):
            scaled = self.fftinterp(subvol[iy,:,:],factor)
            subvol_scaled.append(scaled)
            profs.append(np.mean(scaled,axis=1))


        sy = len(profs[0])
        
        shifts = [0]
        for prof in profs[1:]:
            tx,g = translation1(profs[0],prof)
            shifts.append(tx)

        shifts = -np.array(shifts)
        
        shifts = shifts - np.min(shifts)
        cropped_sy = sy - np.max(shifts)

        aligned = []
        
        for im,x in zip(subvol_scaled,shifts):
            im = im[x:x+cropped_sy,:]
            aligned.append(im)

        aligned = np.array(aligned)
        
        proj = np.mean(aligned,axis=0)

        # fix vertical motion
        
        plt.imshow(proj,cmap='gray')
        plt.show()
        

    def fftinterp(self,im,factor=[2,2]):
        if len(factor)==1:
            factor = [factor,factor]
        if factor[0]==1 and factor[1]==1:
            resized = im
        else:
            isy,isx = im.shape
            yfactor = factor[0]
            xfactor = factor[1]

            osy,osx = int(round(isy*yfactor)),int(round(isx*xfactor))
            while osy%2:
                osy = osy+1
            while osx%2:
                osx = osx+1

            #resized = np.real(np.fft.ifft2(np.fft.fftshift(np.fft.fft2(im)),s=(osy,osx)))

            temp = np.fft.fftshift(np.fft.fft2(im))
            resized = np.zeros((osy,osx),dtype='complex64')
            y1 = (osy-isy)/2
            y2 = y1+isy
            x1 = (osx-isx)/2
            x2 = x1+isx
            resized[y1:y2,x1:x2] = temp
            resized = np.fft.ifftshift(resized)
            resized = np.fft.ifft2(resized)
            resized = np.abs(resized)

            #resized = np.real(np.fft.ifft2(np.fft.fft2(im),s=(osy,osx)))
        return resized
        
        
        
    def blur(self,bscan,kernel_width=5):
        return sp.signal.convolve2d(bscan,np.ones((1,kernel_width)),mode='same')/float(kernel_width)
        
    def click_collector(self,im,message=''):
        global xclicks,yclicks
        xclicks = []
        yclicks = []
        fig = plt.figure()
        
        def onclick(event):
            global clicks
            xnewclick = round(event.xdata)
            ynewclick = round(event.ydata)
            xclicks.append(xnewclick)
            yclicks.append(ynewclick)
            plt.plot(xnewclick,ynewclick,'gs')
            plt.draw()
            
        cid = fig.canvas.mpl_connect('button_press_event',onclick)
        
        plt.imshow(im,cmap='gray')
        plt.autoscale(False)
        
        plt.title(message)
        plt.show()

        return xclicks,yclicks


    def get_label_dict(self):
        out = {}
        try:
            for key in self.h5.get('model/labels').keys():
                out[key] = self.h5.get('model/labels')[key].value
        except Exception as e:
            print e
        return out

    def clear_labels(self):
        self.h5.delete('model/labels')
        
    def plot_profile(self):
        plt.plot(self.profile)
        label_dict = self.get_label_dict()
        for label in label_dict.keys():
            label_z = label_dict[label]
            th = plt.text(label_z,self.profile[label_z],label,ha='center',va='bottom')
        plt.ylim([np.min(self.profile),np.max(self.profile)*1.1])


    def write_volume_labels(self,z_range=2,opening_strel=(1,15),goodness_threshold=0.25,use_fit=True):
        nvol,nslow,ndepth,nfast = self.h5.get(self.data_block).shape
        if use_fit:
            self.logger.info('write_volume_labels: using fit')
            soffset_matrix = np.round(self.h5.get('model/z_offset_fit')[:]).astype(np.float64)
        else:
            offset_matrix = self.h5.get('model/z_offsets')[:].astype(np.float64)
            self.logger.info('write_volume_labels: using offsets (no fit)')
            soffset_matrix = np.zeros_like(offset_matrix)
            for k in range(offset_matrix.shape[0]):
                soffset_matrix[k,:,:] = grey_opening(offset_matrix[k,:,:],opening_strel)
            

        #goodness_matrix = self.h5.get('model/z_offset_goodness')[:]

        label_keys = self.h5.get('model/labels').keys()
        labels = {}
        volume_labels = {}
        projections = {}
        
        for key in label_keys:
            labels[key] = self.h5.get('model/labels')[key].value
            volume_labels[key] = np.zeros((nvol,nslow,nfast))
            projections[key] = np.zeros((nvol,nslow,nfast))
            
        for ivol in range(nvol):
            avol = np.abs(self.h5.get(self.data_block)[ivol,:,:,:])
            self.logger.info('write_volume_labels: Labeling volume %d of %d.'%(ivol+1,nvol))
            self.logger.info('write_volume_labels: Labels: %s.'%','.join(label_keys))

            for islow in range(nslow):
                if (islow+1)%20==0:
                    self.logger.info('write_volume_labels: %d percent done.'%(float(islow+1)/float(nslow)*100))
                for ifast in range(nfast):
                    test = avol[islow,:,ifast]
                    offset = soffset_matrix[ivol,islow,ifast]
                    for key in label_keys:
                        model_z_index = labels[key]
                        volume_labels[key][ivol,islow,ifast] = model_z_index-offset
                        projections[key][ivol,islow,ifast] = np.mean(test[model_z_index-offset-z_range:model_z_index-offset+z_range+1])

        for key in label_keys:
            location = 'model/volume_labels/%s'%key
            plocation = 'projections/%s'%key
            self.h5.put(location,volume_labels[key])
            self.h5.put(plocation,projections[key])
        
                
    def write_axial_alignment(self):
        self.logger.info('write_axial_alignment: Starting.')
        nvol,nslow,ndepth,nfast = self.h5.get(self.data_block).shape
        offset_matrix = self.h5.make('model/z_offsets',(nvol,nslow,nfast),dtype='i2')
        goodness_matrix = self.h5.make('model/z_offset_goodness',(nvol,nslow,nfast),dtype='f8')
        fit_matrix = self.h5.make('model/z_offset_fit',(nvol,nslow,nfast),dtype='f8')
        om,gm,fm = self.align_volumes()
        self.logger.info('write_axial_alignment: Writing axial offsets, goodness, and fit to data store.')
        offset_matrix[...] = om[...]
        goodness_matrix[...] = gm[...]
        fit_matrix[...] = fm[...]
        
    def align_volumes(self):

        # first, we need to make some empty matrices of the right sizes,
        # one layer for each volume, plus the lateral dimensions nslow and nfast
        nvol,nslow,ndepth,nfast = self.h5.get(self.data_block).shape
        offset_matrix = np.zeros((nvol,nslow,nfast))
        goodness_matrix = np.zeros((nvol,nslow,nfast))
        fit_surface_matrix = np.zeros((nvol,nslow,nfast))

        for ivol in range(nvol):
            self.logger.info('align_volumes: working on volume %d of %d.'%(ivol+1,nvol))
            #offset,goodness,fit = self.align_volume(vidx=ivol)
            offset,goodness,fit = self.align_volume_multiscale(vidx=ivol)

            
            offset_matrix[ivol,:,:] = offset
            goodness_matrix[ivol,:,:] = goodness
            fit_surface_matrix[ivol,:,:] = fit

        return offset_matrix,goodness_matrix,fit_surface_matrix


    def align_volume_multiscale(self,vidx=0,rad=2,debug=True):
        self.logger.info('align_volume_multiscale: Starting')
        self.logger.info('align_volume_multiscale: Getting volume from data store.')
        avol = np.abs(self.h5.get(self.data_block)[vidx,:,:,:])
        avol = np.swapaxes(avol,0,1)


        avol = (avol-np.mean(avol,axis=0))/np.std(avol,axis=0)
        

        # NOTE TO SELF: DON'T FUCKING FORGET TO REMOVE THIS, YOU FUCKING IDIOT
        # IF YOU DO, IT WILL RUIN YOUR LIFE FOR A MONTH WHILE YOU DEBUG IT
        # NOTE TO SELF: DON'T FUCKING FORGET TO REMOVE THIS, YOU FUCKING IDIOT
        # IF YOU DO, IT WILL RUIN YOUR LIFE FOR A MONTH WHILE YOU DEBUG IT
        # NOTE TO SELF: DON'T FUCKING FORGET TO REMOVE THIS, YOU FUCKING IDIOT
        # IF YOU DO, IT WILL RUIN YOUR LIFE FOR A MONTH WHILE YOU DEBUG IT
        # NOTE TO SELF: DON'T FUCKING FORGET TO REMOVE THIS, YOU FUCKING IDIOT
        # IF YOU DO, IT WILL RUIN YOUR LIFE FOR A MONTH WHILE YOU DEBUG IT
        # NOTE TO SELF: DON'T FUCKING FORGET TO REMOVE THIS, YOU FUCKING IDIOT
        # IF YOU DO, IT WILL RUIN YOUR LIFE FOR A MONTH WHILE YOU DEBUG IT
        # NOTE TO SELF: DON'T FUCKING FORGET TO REMOVE THIS, YOU FUCKING IDIOT
        # IF YOU DO, IT WILL RUIN YOUR LIFE FOR A MONTH WHILE YOU DEBUG IT
        avol = avol[:,:-2,:]

        

        ndepth,nslow,nfast = avol.shape
        profile = np.zeros(self.profile.shape)
        profile[...] = self.profile[...]
        profile = (profile-np.mean(profile))/np.std(profile)


        
        if len(profile)>ndepth:
            profile = profile[:ndepth]
        if ndepth>len(profile):
            avol = avol[:len(profile),:,:]
            ndepth,nslow,nfast = avol.shape

        gross_profile = np.mean(np.mean(avol,axis=2),axis=1)
        gross_offset,gross_goodness = translation1(gross_profile,profile,debug=False)
        
        offset_submatrix = np.ones((nslow,nfast))*gross_offset
        goodness_submatrix = np.ones((nslow,nfast))*gross_goodness

        # def smooth(v,k):
        #     # an interface to lateral_smooth_3d that handles the
        #     # transposing to minimize confusion
            
        #     return np.transpose(lateral_smooth_3d(np.transpose(avol_bc,(2,0,1)),k),(1,2,0))

        def show_brightest_layer(v):
            # assumes depth is first dimension
            p = np.mean(np.mean(v,axis=2),axis=1)
            mdepth = np.argmax(p)
            layer = v[mdepth,:,:]
            plt.figure()
            plt.imshow(layer)
            plt.colorbar()

        

        # fft the model, for cross-correlation by broadcasting
        f0 = np.fft.fft(profile,axis=0)
        self.logger.info('align_volume_multiscale: model length %d, FFT length %d.'%(len(profile),len(f0)))

        started = False
        sz,sy,sx = avol.shape
        initial_step = float(max(sy,sx))
        final_exp = -np.log(1.49999/initial_step)
        steps = np.round(200*np.exp(-np.linspace(0.0,final_exp,10)))
        
        for k in steps:
            self.logger.info('align_volume_multiscale: smoothing layers in the volume')
            smoothed_avol = lateral_smooth_3d(avol,k)
            # depth is still first dimension

            f1 = np.fft.fft(smoothed_avol,axis=0)
            f1c = f1.conjugate()

            # transpose depth to the final dimension:
            f1_transposed = np.transpose(f1,(1,2,0))
            f1c_transposed = np.transpose(f1c,(1,2,0))
            num_transposed = f0*f1c_transposed
            denom_transposed = np.abs(f0)*np.abs(f1_transposed)
            frac_transposed = num_transposed/denom_transposed

            # now transpose back:
            frac = np.transpose(frac_transposed,(2,0,1))
            
            
            ir = np.abs(np.fft.ifft(frac,axis=0))
            print np.mean(ir),np.max(ir),np.min(ir)

            
            goodness = np.max(ir,axis=0)
            tx = np.argmax(ir,axis=0)

            if not started:
                #goodness_final = goodness.copy()
                #tx_final = tx.copy()
                goodness_final = np.zeros(goodness.shape)
                tx_final = np.zeros(tx.shape)

                goodness_final[...] = goodness[...]
                tx_final[...] = tx[...]
                
                started = True
                better_goodness_mask = np.zeros(goodness.shape)
                small_shifts_mask = np.zeros(tx.shape)
            else:
                better_goodness = (goodness-goodness_final)>-.2
                shift_limit = 10
                small_shifts = np.abs(tx-tx_final)<shift_limit

                better_goodness_mask = np.zeros(goodness.shape)
                better_goodness_mask[np.where(better_goodness)] = 1.0
                small_shifts_mask = np.zeros(tx.shape)
                small_shifts_mask[np.where(small_shifts)] = 1.0


            refinements = better_goodness_mask*small_shifts_mask
            refined_idx = np.where(refinements)
            print refined_idx
            goodness_final[refined_idx] = goodness[refined_idx]
            tx_final[refined_idx] = tx[refined_idx]
            
            if debug:
                plt.clf()
                plt.subplot(2,4,1)
                plt.cla()
                plt.imshow(goodness-goodness_final,cmap='gray')
                plt.colorbar()
                plt.title('$\Delta$ current goodness k=%d'%k)

                plt.subplot(2,4,2)
                plt.cla()
                plt.imshow(tx-tx_final,cmap='gray')
                plt.colorbar()
                plt.title('$\Delta$ current tx k=%d'%k)

                plt.subplot(2,4,3)
                plt.cla()
                plt.imshow(better_goodness_mask,cmap='gray')
                plt.title('better_goodness_mask')
                plt.colorbar()
                
                plt.subplot(2,4,4)
                plt.cla()
                plt.imshow(refinements,cmap='gray')
                plt.title('refinements')
                plt.colorbar()
                
                plt.subplot(2,4,5)
                plt.cla()
                plt.imshow(goodness_final,cmap='gray')
                plt.colorbar()
                plt.title('overall goodness k=%d'%k)

                plt.subplot(2,4,6)
                plt.cla()
                plt.imshow(tx_final,cmap='gray')
                plt.colorbar()
                plt.title('overall tx k=%d'%k)

                plt.subplot(2,4,7)
                plt.cla()
                plt.imshow(small_shifts_mask,cmap='gray')
                plt.colorbar()
                plt.title('small_shifts_mask')
                
                plt.subplot(2,4,8)
                plt.cla()
                plt.imshow(np.mean(smoothed_avol,axis=2),cmap='gray')
                plt.colorbar()

                
                plt.show()

            

    
    def align_volume(self,vidx=0,rad=5):
        self.logger.info('align_volume: Starting')
        self.logger.info('align_volume: Getting volume from data store.')
        avol = np.abs(self.h5.get(self.data_block)[vidx,:,:,:])
        avol = np.swapaxes(avol,0,1)

        

        self.logger.info('align_volume: Smoothing volume with smoothing kernel of size %d pixels.'%rad)
        if rad:
            print time.time()
            avol = lateral_smooth_3d(avol,rad)
            print time.time()
            sys.exit()
            
        ndepth,nslow,nfast = avol.shape
        offset_submatrix = np.zeros((nslow,nfast))
        goodness_submatrix = np.zeros((nslow,nfast))
        profile = self.profile
        
        if len(profile)>ndepth:
            profile = profile[:ndepth]
        if ndepth>len(profile):
            avol = avol[:len(profile),:,:]
            ndepth,nslow,nfast = avol.shape

        x = []
        y = []
        z = []
        w = []
        for islow in range(nslow):
            pct_done = int(round(100*float(islow)/float(nslow)))
            if islow%10==0:
                self.logger.info('align_volume: Aligning A-scans, volume %d is %d percent done.'%(vidx,pct_done))
            for ifast in range(nfast):
                test = avol[:,islow,ifast]
                offset,goodness = translation1(profile,test,debug=False)
                x.append(ifast)
                y.append(islow)
                z.append(offset)
                w.append(goodness)
            
                offset_submatrix[islow,ifast] = offset
                goodness_submatrix[islow,ifast] = goodness

        fitting = True
        if fitting: # revisit this later; may be of use
            ptile = 75
            goodness_threshold=np.percentile(w,ptile)
            self.logger.info('align_volume: Goodness %dth percentile %0.3f used as threshold.'%(ptile,goodness_threshold))
            valid = np.where(w>goodness_threshold)[0]
            x0 = x
            y0 = y
            self.logger.info('align_volume: Using %d of %d points (%d percent) for fit.'%(len(valid),len(w),float(len(valid))/float(len(w))*100))
            x = np.array(x)
            y = np.array(y)
            z = np.array(z)
            w = np.array(w)

            x = x[valid]
            y = y[valid]
            z = z[valid]
            w = w[valid]


            mode='median_filter'
            
            if mode=='spline':
                self.logger.info('Spline fitting surface to A-line axial positions.')
                tck = bisplrep(x,y,z,w=w,xb=0,xe=nfast-1,yb=0,ye=nslow-1)
                self.logger.info('Evaluating spline function at A-line coordinates.')
                fit_surface = bisplev(np.arange(nslow),np.arange(nfast),tck)

            if mode=='polyfit2d':
                self.logger.info('Polynomial fitting surface to A-line axial positions.')
                p = polyfit2d(x,y,z,order=2)
                self.logger.info('Evaluating polynomial function at A-line coordinates.')
                xx,yy = np.meshgrid(np.arange(nfast),np.arange(nslow))
                fit_surface = polyval2d(xx,yy,p)

            if mode=='median_filter':
            # This is a dumb way to fit. Use the goodness matrix, dummy, perhaps with 2D splines!
                self.logger.info('Median filtering to create a smoothed offset surface.')
                slow_height = np.median(offset_submatrix,axis=1)
                plt.figure()
                plt.plot(slow_height)
                plt.figure()
                debias = (offset_submatrix.T-slow_height).T

                
                plt.figure()
                plt.imshow(debias)
                plt.colorbar()
                
                fit_surface_1 = median_filter(offset_submatrix,(3,3))
                plt.figure()
                plt.imshow(fit_surface_1)
                plt.title('straight fit')
                plt.colorbar()
                
                fit_surface_2 = (median_filter(debias,(3,3)).T+slow_height).T
                plt.figure()
                plt.imshow(fit_surface_2)
                plt.title('debiased fit')
                plt.colorbar()
                plt.show()
                sys.exit()

            if mode=='interp2d':
                self.logger.info('Using interp2d to create a smoothed offset surface.')
                interpolation_function = interp2d(x,y,z)
                fit_surface = interpolation_function(x0,y0)
                print fit_surface
                print fit_surface.shape

            
    
            # print fit_surface
            # print fit_surface.shape
            if True:
                clim = np.min(offset_submatrix),np.max(offset_submatrix)
                plt.figure()
                plt.imshow(offset_submatrix,interpolation='none',clim=clim)
                plt.colorbar()
                plt.figure()
                plt.imshow(fit_surface,interpolation='none',clim=clim)
                plt.colorbar()
                plt.figure()
                plt.imshow(offset_submatrix-fit_surface,interpolation='none')
                plt.colorbar()

                plt.show()


            #goodness_used = np.zeros(goodness_submatrix.shape)
            #goodness_used[np.where(goodness_submatrix>goodness_threshold)] = 1.0
        else:
            fit_surface = offset_submatrix
            
        return offset_submatrix,goodness_submatrix,fit_surface
                
                
    def cleanup_modeldb(self):
        modeldb = H5(ocfg.model_database)
        for key in modeldb.h5.keys():
            if len(modeldb.h5[key]['labels'].keys())==0:
                del modeldb.h5[key]
        modeldb.close()
        
    def find_matching_labels(self,goodness_threshold=0.6,rad=3):
        modeldb = H5(ocfg.model_database)
        print modeldb.h5.keys()
        guess_dicts = [{}]
        goodnesses = [-np.inf]
        for key in modeldb.h5.keys():
            guess_dict = {}
            test = modeldb.get(key)['profile'][:]
            if len(test)>len(self.profile):
                test_profile = test[:len(self.profile)]
            elif len(test)<len(self.profile):
                test_profile = np.zeros(len(self.profile))
                test_profile[:len(test)] = test
            else:
                test_profile = test
            offset,goodness = translation1(test_profile,self.profile,debug=False)
            self.logger.info('Comparing with model labels from %s: goodness %0.3f.'%(key,goodness))
            if goodness>goodness_threshold:
                test_labels = modeldb.get(key)['labels'].keys()
                if len(test_labels)>0:
                    for test_label in test_labels:
                        search_position = modeldb.get(key)['labels'][test_label].value - offset
                        temp = np.zeros(self.profile.shape)
                        temp[search_position-rad:search_position+rad+1] = self.profile[search_position-rad:search_position+rad+1]
                        guess_dict[test_label] = np.argmax(temp)
                    guess_dicts.append(guess_dict)
                    goodnesses.append(goodness)
        
        return guess_dicts[np.argmax(goodnesses)]

    
    def click_label(self,smoothing=1):
        if smoothing>1:
            working_profile = sp.signal.convolve(self.profile,np.ones((smoothing)),mode='same')/float(smoothing)
        else:
            working_profile = self.profile
        
        # find peaks and troughs:
        gthresh = 1.0/smoothing

        nslow = self.h5.h5[self.data_block].shape[1]
        

        peaks = np.sort(find_peaks(working_profile,gradient_threshold=gthresh))
        
        idx = 0
        z = np.arange(len(working_profile))
        done = False or not len(peaks)


        fig = plt.figure(figsize=(22,12))
        for key in plt.rcParams.keys():
            if key[:6]=='keymap':
                plt.rcParams[key] = ''
        
        global current_x,current_label,label_dict

        # try to get the label dictionary from the current dataset
        # if none exists, mine the model database for a match
        
        label_dict = self.get_label_dict()
        if len(label_dict)==0:
            label_dict = self.find_matching_labels()
            
        current_x = 0
        current_label = ''

        l1 = .05
        l2 = .55
        fw = .9
        hw = .4
        b1 = .55
        b2 = .05
        fh = .9
        hh = .4
        
        global bscanindex 
        bscanindex = 0

        def plot_at(x):
            global current_label,bscanindex
            
            if x in label_dict.values():
                existing_label = [key for key, value in label_dict.items() if value == x][0]
            else:
                existing_label = ''
            
            plt.axes([l1,b2,hw,hh])
            plt.cla()
            bscan = np.abs(self.h5.h5[self.data_block][0,bscanindex,:,:])
            #bscan = shear(bscan,1)
            try:
                test = np.mean(bscan[:,-20:],axis=1)
            except:
                test = np.mean(bscan,axis=1)
            offset,goodness = translation1(test,working_profile,xlims=10,equalize=True)
            
            cmin = np.median(bscan)
            cmax = np.percentile(bscan,99.95) # saturate 0.05% of pixels
            plt.imshow(bscan,interpolation='none',clim=(cmin,cmax),cmap='gray')
            for label in label_dict.keys():
                print label,label_dict[label]
                label_z = z[label_dict[label]]
                th = plt.text(bscan.shape[1],label_z-offset,label,ha='left',va='center',fontsize=8)
            try:
                plt.ylim((np.max(label_dict.values())+10,np.min(label_dict.values())-10))
            except:
                pass
                
            plt.axes([l1,b1,fw,hh])
            plt.cla()
            plt.plot(z,working_profile)
            plt.plot(z[x],working_profile[x],'ks')
            valid = np.where(working_profile)[0]
            plt.xlim((valid[0],valid[-1]))
            plt.autoscale(False)
            for label in label_dict.keys():
                label_z = z[label_dict[label]]
                th = plt.text(label_z,working_profile[label_z],label,ha='center',va='bottom')
                
            plt.axes([l2,b2,hw,hh])
            plt.cla()
            plt.plot(z,working_profile)
            plt.plot(z[x],working_profile[x],'ks')
            plt.xlim((z[x]-10,z[x]+10))
            z1 = max(0,x-10)
            z2 = min(len(z),x+10)
            ymin = np.min(working_profile[z1:z2])
            ymax = np.max(working_profile[z1:z2])*1.25
            plt.text(z[x],working_profile[x],existing_label,ha='center',va='bottom')
            plt.ylim((ymin,ymax))
            plt.title(current_label)
            plt.draw()


        def onclick(event):
            global current_x
            # print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
            #     event.button, event.x, event.y, event.xdata, event.ydata)
            current_x = round(event.xdata)
            plot_at(current_x)

        def onpress(event):
            #print event.key
            global current_x,current_label,label_dict,bscanindex
            if event.key=='right':
                current_x = (current_x + 1)%len(working_profile)
            elif event.key=='ctrl+right':
                try:
                    current_x = np.min(peaks[np.where(peaks>current_x)[0]])
                except Exception as e:
                    current_x = len(working_profile)-1
            elif event.key=='left':
                current_x = (current_x - 1)%len(working_profile)
            elif event.key=='ctrl+left':
                try:
                    current_x = np.max(peaks[np.where(peaks<current_x)[0]])
                except Exception as e:
                    current_x = 0
            elif event.key=='shift+ctrl+right':
                try:
                    current_x = peaks[np.where(peaks>current_x)[0]][5]
                except Exception as e:
                    current_x = len(working_profile)-1
            elif event.key=='shift+ctrl+left':
                try:
                    current_x = peaks[np.where(peaks<current_x)[0]][-5]
                except Exception as e:
                    current_x = 0
            elif event.key=='shift':
                pass
            elif event.key=='/':
                pass
            elif event.key in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ':
                current_label = current_label + event.key.upper()
            elif event.key=='backspace':
                current_label = current_label[:-1]
            elif event.key=='ctrl+delete':
                label_dict = {}
            elif event.key=='delete':
                for key in label_dict.keys():
                    if label_dict[key]==current_x:
                        label_dict.pop(key)
            elif event.key=='enter':
                label_dict[current_label] = current_x
                print label_dict
                current_label = ''
            elif event.key=='pageup':
                bscanindex = (bscanindex + 1)%nslow
            elif event.key=='pagedown':
                bscanindex = (bscanindex - 1)%nslow
                
            plot_at(current_x)
            

        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        pid = fig.canvas.mpl_connect('key_press_event', onpress)
        
        plot_at(current_x)
        
        plt.show()

        self.h5.require_group('model')
        self.h5.require_group('model/labels')
        for key in label_dict.keys():
            self.h5.put('model/labels/%s'%key,label_dict[key])
            
        mdb = H5(ocfg.model_database)
        did = self.h5.get('IDs/dataset_id').value
        # did is the primary key for the model, but we'll also save eccentricity
        did_key = '%d'%did
        mdb.require_group(did_key)

        si = self.h5.get('eccentricity/superior_inferior').value
        nt = self.h5.get('eccentricity/nasal_temporal').value
        
        radial_distance = np.sqrt(si**2+nt**2)
        
        mdb.put('%s/superior_inferior'%did_key,si)
        mdb.put('%s/nasal_temporal'%did_key,nt)
        mdb.put('%s/radial_distance'%did_key,radial_distance)
        mdb.put('%s/profile'%did_key,self.profile)

        mdb.require_group('%s/labels'%did_key)
        for key in label_dict.keys():
            mdb.put('%s/labels/%s'%(did_key,key),label_dict[key])

        mdb.close()
        


def test():
    h5 = H5('./oct_test_volume/oct_test_volume_2T.hdf5')
    m = Model(h5,True)
    #m.make_model()
    #m.clear_labels()
    #m.click_label()
    m.write_axial_alignment()
    m.write_volume_labels()
    
if __name__=='__main__':
    test()
