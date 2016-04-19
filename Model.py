import sys,os
import logging
import numpy as np
import scipy as sp
from scipy.ndimage.filters import generic_filter
from scipy.interpolate import bisplrep,bisplev
from matplotlib import pyplot as plt
from utils import translation,translation1,autotrim_bscan,find_peaks,shear,Clock,lateral_smooth_3d
from octopod.Misc import H5
import logging
import octopod_config as ocfg
logging.basicConfig(level=logging.DEBUG)


if False:
    x = np.random.rand(100,100)
    a = x[5:,:]
    b = x[:-5,:]
    print translation(a,b)

class Model:
    """A model of gross retinal reflectance, used for segmenting
    and labeling A-scans."""

    def __init__(self,h5,debug=False):
        """Initialize model. May pass an h5py.File object or .hdf5 filename.
        The file or object must contain a 'processed_data' dataset containing
        at least one volume."""
        self.logger = logging.getLogger(__name__)
        if type(h5)==str:
            self.h5 = H5(h5)
            self.logger.info('Opening file %s.'%self.h5)
        else:
            self.h5 = h5
        try:
            self.profile = self.h5.get('model/profile')[:]
        except Exception as e:
            self.logger.info('Model does not exist in h5 file.')
            self.profile = self.make_model(debug=debug)

    def blur(self,bscan,kernel_width=5):
        return sp.signal.convolve2d(bscan,np.ones((1,kernel_width)),mode='same')/float(kernel_width)
        
    def make_model(self,vidx=0,debug=False):
        self.logger.info('Making model...')
        avol = np.abs(self.h5.get('processed_data')[vidx,:,:,:])
        nSlow,nFast,nDepth = avol.shape
        
        if False:
            # make a test volume with very clean translations.
            avol = np.ones(avol.shape)
            zpos = 100
            rad = 3
            zpos_vec = []
            for iSlow in range(nSlow):
                avol[iSlow,zpos-rad:zpos+rad,:] = avol[iSlow,zpos-rad:zpos+rad,:] + 100.0
                zpos = np.abs(np.round(zpos + np.random.randn()))
                zpos_vec.append(zpos)

                
        # crop to avoid edge effects in the registration:
        if avol.shape[2]>10:
            avol = avol[:,:,5:-5]

        tx_vec = [0.0]
        ty_vec = [0.0]
        g_vec = [1.0]

        template = None

        # compute an intensity threshold for the first template image,
        # to avoid using 

        thresh = np.mean(avol)*.75

        for iSlow in range(1,nSlow):
            last = self.blur(avol[iSlow-1,:,:])

            bright_enough = np.mean(last)>=thresh
            
            if debug:
                self.logger.debug('Frame %d of %d. Sufficiently bright: %d'%(iSlow,nSlow-1,bright_enough))
            
            if template is None:
                if bright_enough:
                    template = last
                    counter = np.ones(template.shape)
                    sy,sx = template.shape
                else:
                    continue
            
            target = self.blur(avol[iSlow,:,:])
            
            if False:
                plt.subplot(1,2,1)
                plt.cla()
                plt.imshow(template,interpolation='none',aspect='auto')
                #plt.colorbar()
                plt.subplot(1,2,2)
                plt.cla()
                plt.imshow(target,interpolation='none',aspect='auto')
                #plt.colorbar()
                plt.pause(.1)
            
            tx,ty,g = translation(target,template,debug=False)
            if debug:
                self.logger.debug('x:%d; y:%d; goodness:%0.3f'%(tx,ty,g))
                
            # ty is the amount to shift target to align with template
            # let's say template's profile is [0 1 0 0] and target's is
            # [0 0 1 0]. ty = -1 in this case.
            if ty<0:
                put0 = -ty
                put1 = sy
                get0 = 0
                get1 = ty
            elif ty>0:
                put0 = 0
                put1 = -ty
                get0 = ty
                get1 = sy
            else:
                get0 = 0
                get1 = sy
                put0 = 0
                put1 = sy

            template[put0:put1,:] = template[put0:put1,:] + target[get0:get1,:]
            counter[put0:put1,:] = counter[put0:put1] + 1.0

            if debug:
                plt.cla()
                plt.plot(np.mean(template/counter,axis=1))
                plt.pause(.1)


        template = template/counter
        template = shear(template,2)

        model_profile = np.mean(template,axis=1)
        
        self.write_profile(model_profile)
        return model_profile

    def write_profile(self,model_profile):
        self.h5.require_group('model')
        self.h5.put('model/profile',model_profile)

    def click_crop(self):
        global clicks
        clicks = []
        
        fig = plt.figure(figsize=(10,6))
        
        def onclick(event):
            global clicks
            newclick = round(event.xdata)
            clicks.append(newclick)
            plt.axvline(newclick)
            plt.draw()
            
        cid = fig.canvas.mpl_connect('button_press_event',onclick)
        
        plt.plot(self.profile)
        plt.show()

        if len(clicks)>=2:

            z1 = np.min(clicks)
            z2 = np.max(clicks)
            self.profile[:z1] = 0.0
            self.profile[z2+1:] = 0.0
            self.write_profile(self.profile)
            


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


    def write_volume_labels(self,z_tolerance=2,medfilt_kernel=(1,5,15),goodness_threshold=0.25):
        nvol,nslow,ndepth,nfast = self.h5.get('processed_data').shape

        offset_matrix = self.h5.get('model/z_offsets')[:].astype(np.float64)
        offset_matrix = np.round(self.h5.get('model/z_offset_fit')[:]).astype(np.float64)

        goodness_matrix = self.h5.get('model/z_offset_goodness')[:]

        foffset_matrix = np.zeros(offset_matrix.shape)
        foffset_matrix[:,:,:] = offset_matrix[:,:,:]
        foffset_matrix[np.where(goodness_matrix<goodness_threshold)] = np.nan
        foffset_matrix = generic_filter(foffset_matrix,np.nanmedian,medfilt_kernel,mode='nearest')
        
        # now, filter out clear outliers; replace with global mode
        mode = sp.stats.mode(offset_matrix.ravel())[0][0]
        
        mask = np.zeros(foffset_matrix.shape)
        lower_threshold = np.mean(offset_matrix)-1.5*np.std(offset_matrix)
        upper_threshold = np.mean(offset_matrix)+1.5*np.std(offset_matrix)
        cond = np.logical_and(foffset_matrix>lower_threshold,foffset_matrix<upper_threshold)
        mask[np.where(cond)] = 1
        self.logger.info('Fraction of pixels in offset matrix deemed valid: %0.4f.'%(np.sum(mask)/np.prod(mask.shape)))
        self.logger.info('Setting invalid pixels to offset mode: %d'%mode)
        foffset_matrix[np.where(1-mask)] = mode
        label_keys = self.h5.get('model/labels').keys()
        labels = {}
        volume_labels = {}
        projections = {}
        
        for key in label_keys:
            labels[key] = self.h5.get('model/labels')[key].value
            volume_labels[key] = np.zeros((nvol,nslow,nfast))
            projections[key] = np.zeros((nvol,nslow,nfast))
            
        for ivol in range(nvol):
            avol = np.abs(self.h5.get('processed_data')[ivol,:,:,:])
            # plt.subplot(121)
            # plt.imshow(offset_matrix[ivol,:,:])
            # plt.colorbar()
            # plt.subplot(122)
            # plt.imshow(foffset_matrix[ivol,:,:])
            # plt.colorbar()
            # plt.show()
            self.logger.info('Labeling volume %d of %d.'%(ivol+1,nvol))

            for islow in range(nslow):
                if (islow+1)%20==0:
                    self.logger.info('%d percent done.'%(float(islow+1)/float(nslow)*100))
                for ifast in range(nfast):
                    test = avol[islow,:,ifast]
                    offset = offset_matrix[ivol,islow,ifast]
                    foffset = foffset_matrix[ivol,islow,ifast]
                    for key in label_keys:
                        model_z_index = labels[key]
                        temp = np.zeros(test.shape)
                        temp[...] = test[...]
                        z1,z2 = model_z_index-offset-z_tolerance,model_z_index-offset+z_tolerance
                        temp[:z1] = -np.inf
                        temp[z2+1:] = -np.inf
                        
                        # plt.subplot(311)
                        # plt.cla()
                        # self.plot_profile()
                        # plt.subplot(312)
                        # plt.cla()
                        # plt.plot(test)
                        # plt.axvspan(z1,z2,alpha=0.3)
                        # plt.title(key)
                        # print offset,z1,z2
                        # plt.subplot(313)
                        # plt.cla()
                        # plt.plot(temp)
                        # plt.pause(1)
                        max_idx = np.argmax(temp)
                        volume_labels[key][ivol,islow,ifast] = max_idx
                        projections[key][ivol,islow,ifast] = np.mean(test[max_idx-z_tolerance:max_idx+z_tolerance+1])
            print

        for key in label_keys:
            location = 'model/volume_labels/%s'%key
            plocation = 'projections/%s'%key
            self.h5.put(location,volume_labels[key])
            self.h5.put(plocation,projections[key])
        
                
    def write_axial_alignment(self):
        nvol,nslow,ndepth,nfast = self.h5.get('processed_data').shape
        offset_matrix = self.h5.make('model/z_offsets',(nvol,nslow,nfast),dtype='i2')
        goodness_matrix = self.h5.make('model/z_offset_goodness',(nvol,nslow,nfast),dtype='f8')
        fit_matrix = self.h5.make('model/z_offset_fit',(nvol,nslow,nfast),dtype='f8')
        om,gm,fm = self.align_volumes()
        offset_matrix[...] = om[...]
        goodness_matrix[...] = gm[...]
        fit_matrix[...] = fm[...]
        
    def align_volumes(self):
        nvol,nslow,ndepth,nfast = self.h5.get('processed_data').shape
        offset_matrix = np.zeros((nvol,nslow,nfast))
        goodness_matrix = np.zeros((nvol,nslow,nfast))
        fit_surface_matrix = np.zeros((nvol,nslow,nfast))
        
        for ivol in range(nvol):
            offset,goodness,fit = self.align_volume(vidx=ivol)
            offset_matrix[ivol,:,:] = offset
            goodness_matrix[ivol,:,:] = goodness
            fit_surface_matrix[ivol,:,:] = fit

        return offset_matrix,goodness_matrix,fit_surface_matrix
        
    def align_volume(self,vidx=0,rad=5,goodness_threshold=0.0):
        avol = np.abs(self.h5.get('processed_data')[vidx,:,:,:])
        avol = np.swapaxes(avol,0,1)

        if rad:
            avol = lateral_smooth_3d(avol,rad)
        
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
            self.logger.info('B-scan %d in volume %d.'%(islow,vidx))
            for ifast in range(nfast):
                #self.logger.info('A-scan %d.'%ifast)
                test = avol[:,islow,ifast]
                offset,goodness = translation1(profile,test,debug=False)
                if goodness>goodness_threshold:
                    x.append(ifast)
                    y.append(islow)
                    z.append(offset)
                    w.append(goodness)
                
                offset_submatrix[islow,ifast] = offset
                goodness_submatrix[islow,ifast] = goodness

        tck = bisplrep(x,y,z,w=w,xb=0,xe=nfast-1,yb=0,ye=nslow-1)
        fit_surface = bisplev(np.arange(nfast),np.arange(nslow),tck)

        goodness_used = np.zeros(goodness_submatrix.shape)
        goodness_used[np.where(goodness_submatrix>goodness_threshold)] = 1.0
        
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

        nslow = self.h5.h5['processed_data'].shape[1]
        

        peaks = np.sort(find_peaks(working_profile,gradient_threshold=gthresh))
        
        idx = 0
        z = np.arange(len(working_profile))
        done = False or not len(peaks)


        fig = plt.figure(figsize=(22,16))
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
            bscan = np.abs(self.h5.h5['processed_data'][0,bscanindex,:,:])
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
