import numpy as np
from matplotlib import pyplot as plt
import sys,os
from octopod import H5,utils
import glob
from scipy.ndimage import zoom
from scipy.interpolate import griddata
from scipy.signal import fftconvolve

class Series:

    def __init__(self,series_filename,reference_frame=None):
        self.h5 = H5(series_filename)
        self.n_frames = 0
        self.tag_template = '%s/%03d'
        self.series_filename = series_filename
        self.working_directory = os.path.split(series_filename)[0]
        if not reference_frame is None:
            self.reference = reference_frame
            self.h5.put('/reference_frame',reference_frame)

    def set_reference_frame(self,reference_frame):
        self.reference = reference_frame
        self.h5.put('/reference_frame',reference_frame)



    def correct_reference(self,kernel_size=10):
        try:
            si = self.h5['sum_image'][:,:]
            ai = self.h5['average_image'][:,:]
            ci = self.h5['counter_image'][:,:]

            rx1 = self.h5['reference_coordinates/x1'].value
            rx2 = self.h5['reference_coordinates/x2'].value
            ry1 = self.h5['reference_coordinates/y1'].value
            ry2 = self.h5['reference_coordinates/y2'].value
        except Exception as e:
            sys.exit('Has this Series been rendered? If not, do that first.')

        sy,sx = si.shape
        
        def hcentroid(im):
            hidx = np.arange(im.shape[1])
            return np.sum(im*hidx,axis=1)/np.sum(im,axis=1)
    
        hc = -hcentroid(ci[int(ry1):int(ry2),int(rx1):int(rx2)])
        hc = (hc - np.mean(hc))


        rsy = len(hc)
        cut = rsy//20
        hc[:cut] = hc.mean()
        hc[-cut:] = hc.mean()
        
        def velocity_to_position(vec):
            y = np.cumsum(vec)
            x = np.arange(len(vec))
            fit = np.polyfit(x,y,1)
            yfit = np.polyval(fit,x)
            return (y - yfit)/float(len(vec))


        vprof = np.max(ci,axis=1)[ry1:ry2]

        plt.plot(vprof)
        plt.show()
        sys.exit()
        
        vc = velocity_to_position(vprof)
        vc[:cut] = vc.mean()
        vc[-cut:] = vc.mean()

        comp = hc + vc*1j
        comp = fftconvolve(comp,np.ones((cut)),mode='same')/float(cut)

        hc = np.real(comp)
        vc = np.imag(comp)
        
        reference = self.h5['reference_frame'][:,:]
        original_sy,original_sx = reference.shape
        # downsample the movement estimates:
        for fk in self.h5['frames'].keys():
            for idx in self.h5['frames'][fk].keys():
                oversample_factor = self.h5['frames'][fk][idx]['oversample_factor'].value
                break

        hc = np.reshape(hc,(original_sy,int(oversample_factor))).mean(axis=1)
        vc = np.reshape(vc,(original_sy,int(oversample_factor))).mean(axis=1)

        to_XX,to_YY = np.meshgrid(np.arange(original_sx),np.arange(original_sy))

        from_XX = (to_XX.T + hc).T
        from_YY = (to_YY.T + vc).T

        to_XX = to_XX.ravel()
        to_YY = to_YY.ravel()
        from_XX = from_XX.ravel()
        from_YY = from_YY.ravel()

        # use scipy.interpolate.griddata parlance for clarity
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html
        points = np.vstack((from_XX,from_YY)).T
        values = reference.ravel()

        corrected_reference = griddata(points,values,(to_XX.ravel(),to_YY.ravel()),method='cubic')
        corrected_reference = np.reshape(corrected_reference,reference.shape)

        self.h5.put('/corrected/reference_frame',corrected_reference)
        self.h5.put('/corrected/horizontal_correction',hc)
        self.h5.put('/corrected/vertical_correction',vc)
        
        
    def add(self,filename,vidx,layer_names=None,overwrite=False,oversample_factor=3,strip_width=3.0,do_plot=False):
        
        print 'Adding %s, volume %d.'%(filename,vidx)
        
        target_tag = self.tag_template%(os.path.split(filename)[1],vidx)

        if self.h5.has('/frames/%s'%target_tag) and not overwrite:
            print 'Series already has entry for %s.'%target_tag
            return

        target = self.get_image(filename,vidx,layer_names)
        reference = self.reference
        y,x,g = utils.strip_register(target,reference,oversample_factor,strip_width,do_plot=do_plot)
        
        self.h5.put('/frames/%s/x_shifts'%target_tag,x)
        self.h5.put('/frames/%s/y_shifts'%target_tag,y)
        self.h5.put('/frames/%s/goodnesses'%target_tag,g)
        self.h5.put('/frames/%s/reference'%target_tag,[0])
        self.h5.put('/frames/%s/oversample_factor'%target_tag,oversample_factor)

    def is_registered(self):
        counter = 0
        try:
            for filename in self.h5['/frames'].keys():
                for fileindex in self.h5['/frames'][filename].keys():
                    counter = counter + 1
        except Exception as e:
            print e
        return counter
                
        
    def get_image(self,filename_stub,vidx,layer_names):
        filename = os.path.join(self.working_directory,filename_stub)
        target_h5 = H5(filename)

        if layer_names is None:
            # if the layer_names list is missing, use the first layer as a default
            # this seems like okay behavior since most times there's only one projection
            # anyway
            layer_names = [target_h5['projections'].keys()[0]]

        if len(layer_names)>1:
            test = target_h5['projections'][layer_names[0]][vidx,:,:]
            n_slow,n_fast = test.shape
            stack = np.zeros((len(layer_names),n_slow,n_fast))
            for idx,layer_name in enumerate(layer_names):
                stack[idx,:,:] = target_h5['projections'][layer_name][vidx,:,:]
            out = np.mean(stack,axis=0)
            del stack
        else:
            out = target_h5['projections'][layer_names[0]][vidx,:,:]    
        return out

    def get_n_frames(self):
        count = 0
        try:
            filenames = self.h5['/frames'].keys()
            for filename in filenames:
                count = count + len(self.h5['/frames/%s'%filename].keys())
        except:
            pass
        return count


    def is_rendered(self):
        out = True
        try:
            test = self.h5['/sum_image']
        except:
            out = False

        return out

            
    def render(self,layer_names=None,goodness_threshold=0.0,correlation_threshold=0.0,overwrite=False,oversample_factor=3,do_plot=False,frames_to_save=[],corrected=False):

        if len(frames_to_save):
            frames_directory = self.series_filename.replace('.hdf5','')+'_saved_frames'
            if not os.path.exists(frames_directory):
                os.makedirs(frames_directory)

        files = self.h5['frames'].keys()
        
        sign = -1
            
        # first, find the minimum and maximum x and y shifts
        xmin = np.inf
        xmax = -np.inf
        ymin = np.inf
        ymax = -np.inf

        for filename in files:
            keys = self.h5['/frames/%s'%filename].keys()
            for k in keys:
                test = self.get_image(filename,0,None)
                n_slow,n_fast = test.shape

                if corrected:
                    hc = self.h5['/corrected/horizontal_correction'][:]
                    vc = self.h5['/corrected/vertical_correction'][:]
                else:
                    hc = np.zeros(n_slow)
                    vc = np.zeros(n_slow)
                
                goodnesses = self.h5['/frames/%s/%s/goodnesses'%(filename,k)][:]
                xshifts = sign*self.h5['/frames/%s/%s/x_shifts'%(filename,k)][:]
                yshifts = sign*self.h5['/frames/%s/%s/y_shifts'%(filename,k)][:]

                xshifts = np.squeeze(xshifts)+hc
                yshifts = np.squeeze(yshifts)+vc

                xshifts,yshifts,goodnesses = self.filter_registration(xshifts,yshifts,goodnesses)

                yshifts = yshifts + np.arange(n_slow)

                assert len(yshifts)==n_slow
                
                valid = np.where(goodnesses>=goodness_threshold)[0]

                if len(valid):
                    xshifts = xshifts[valid]
                    yshifts = yshifts[valid]

                    newxmin = np.min(xshifts)
                    newxmax = np.max(xshifts)
                    newymin = np.min(yshifts)
                    newymax = np.max(yshifts)

                    xmin = min(xmin,newxmin)
                    xmax = max(xmax,newxmax)
                    ymin = min(ymin,newymin)
                    ymax = max(ymax,newymax)


        xmin = np.round(xmin*oversample_factor)
        xmax = np.round(xmax*oversample_factor)
        dx = xmax-xmin
        xoffset = xmin
        
        width = n_fast*oversample_factor + dx

        ymin = np.round(ymin*oversample_factor)
        ymax = np.round(ymax*oversample_factor)
        dy = ymax-ymin
        yoffset = ymin
        height = oversample_factor + dy

        test_oversampled = zoom(test,oversample_factor)
        sy_oversampled,sx_oversampled = test_oversampled.shape

        erx1 = -xoffset
        ery1 = -yoffset
        erx2 = erx1+sx_oversampled
        ery2 = ery1+sy_oversampled
        
        embedded_reference = np.zeros((height,width))
        ref_oversampled = zoom(self.reference,oversample_factor)
        embedded_reference[ery1:ery2,erx1:erx2] = ref_oversampled

        sum_image = np.zeros((height,width))
        counter_image = np.zeros((height,width))
        correlation_image = np.zeros((height,width))
        
        x1 = round(sign*xoffset)
        x2 = x1+sx_oversampled
        y1 = round(sign*yoffset)
        y2 = y1+sy_oversampled
        fig = plt.figure()
        all_corr_coefs = []
        #ref_clim = np.percentile(self.reference,(1,99.5))


        frame_index = 0
        for filename in files:

            keys = self.h5['/frames/%s'%filename].keys()
            
            for k in keys:
                current_sum_image = np.zeros((height,width))
                current_counter_image = np.zeros((height,width))
                goodnesses = self.h5['/frames/%s/%s/goodnesses'%(filename,k)][:]
                xshifts = sign*self.h5['/frames/%s/%s/x_shifts'%(filename,k)][:]
                yshifts = sign*self.h5['/frames/%s/%s/y_shifts'%(filename,k)][:]
                xshifts = np.squeeze(xshifts)+hc
                yshifts = np.squeeze(yshifts)+vc

                xshifts,yshifts,goodnesses = self.filter_registration(xshifts,yshifts,goodnesses)

                valid = np.where(goodnesses>=goodness_threshold)[0]

                im = self.get_image(filename,int(k),layer_names)

                if (not any(xshifts)) and (not any(yshifts)) and (not corrected):
                    block = zoom(im,oversample_factor)
                    bsy,bsx = block.shape
                    x1 = -xoffset
                    y1 = -yoffset
                    x2 = x1+bsx
                    y2 = y1+bsy
                    sum_image[y1:y2,x1:x2] = sum_image[y1:y2,x1:x2] + block
                    counter_image[y1:y2,x1:x2] = counter_image[y1:y2,x1:x2] + 1.0
                    correlation_image[y1:y2,x1:x2] = correlation_image[y1:y2,x1:x2] + 1.0
                    correlation_vector = np.ones(goodnesses.shape)
                
                    current_sum_image[y1:y2,x1:x2] = current_sum_image[y1:y2,x1:x2] + block
                    current_counter_image[y1:y2,x1:x2] = current_counter_image[y1:y2,x1:x2] + 1.0
                    self.h5.put('/reference_coordinates/x1',x1)
                    self.h5.put('/reference_coordinates/x2',x2)
                    self.h5.put('/reference_coordinates/y1',y1)
                    self.h5.put('/reference_coordinates/y2',y2)

                else:
                    correlation_vector = np.ones(goodnesses.shape)*np.nan
                    if len(valid):
                        for v in valid:
                            line = im[v]
                            line = np.expand_dims(line,0)
                            block = zoom(line,oversample_factor)
                            bsy,bsx = block.shape
                            x1 = np.round(xshifts[v]*oversample_factor-xoffset)
                            y1 = v*oversample_factor+np.round(yshifts[v]*oversample_factor-yoffset)
                            x2 = x1+bsx
                            y2 = y1+bsy

                            

                            # im1 = sum_image[y1:y2,x1:x2].ravel()
                            # if np.min(im1)==np.max(im1) and False:
                            #     corr = 1.0
                            # else:
                            #     corr_valid = np.where(im1)[0]
                            #     a,b = im1[corr_valid],block.ravel()[corr_valid]
                            #     corr = np.corrcoef(a,b)[1,0]
                            #     plt.figure()
                            #     plt.subplot(1,2,1)
                            #     plt.plot(a)
                            #     plt.subplot(1,2,2)
                            #     plt.plot(b)
                            #     plt.show()
                            # all_corr_coefs.append(corr)

                            correlation_image[y1:y2,x1:x2] = correlation_image[y1:y2,x1:x2] + corr
                            plt.figure(2)
                            plt.clf()
                            plt.imshow(correlation_image)
                            plt.colorbar()
                            plt.pause(.1)
                            if corr>correlation_threshold:
                                sum_image[y1:y2,x1:x2] = sum_image[y1:y2,x1:x2] + block
                                counter_image[y1:y2,x1:x2] = counter_image[y1:y2,x1:x2] + 1.0
                                current_sum_image[y1:y2,x1:x2] = current_sum_image[y1:y2,x1:x2] + block
                                current_counter_image[y1:y2,x1:x2] = current_counter_image[y1:y2,x1:x2] + 1.0

                if do_plot:
                    temp = counter_image.copy()
                    temp[np.where(temp==0)] = 1.0
                    av = sum_image/temp


                    plt.clf()
                    plt.subplot(1,3,1)
                    plt.cla()
                    clim = (np.median(av)-2*np.std(av),np.median(av)+2*np.std(av))
                    plt.imshow(av,cmap='gray',clim=clim,interpolation='none')

                    plt.subplot(1,3,2)
                    plt.cla()
                    plt.imshow(counter_image)
                    plt.colorbar()

                    plt.subplot(1,3,3)
                    plt.cla()
                    try:
                        plt.hist(all_corr_coefs,100)
                    except:
                        pass
                    #plt.colorbar()
                    plt.pause(.0000000001)


                if frame_index in frames_to_save:
                    if corrected:
                        outfn = os.path.join(frames_directory,'corrected_frame_%03d.npy'%frame_index)
                    else:
                        outfn = os.path.join(frames_directory,'frame_%03d.npy'%frame_index)
                    current_counter_image[np.where(current_counter_image==0)] = 1.0
                    out = current_sum_image/current_counter_image
                    print 'Saving frame to %s.'%outfn
                    np.save(outfn,out)
                frame_index = frame_index + 1
            
        temp = counter_image.copy()
        temp[np.where(temp==0)] = 1.0
        av = sum_image/temp

        if corrected:
            corrstr = '/corrected'
        else:
            corrstr = ''
            
        self.h5.put('%s/counter_image'%corrstr,counter_image)
        self.h5.put('%s/average_image'%corrstr,av)
        self.h5.put('%s/sum_image'%corrstr,sum_image)
        

        if do_plot:
            plt.close()

            plt.subplot(1,2,1)
            plt.imshow(av,cmap='gray',clim=clim,interpolation='none')

            plt.subplot(1,2,2)
            plt.imshow(counter_image)
            plt.colorbar()

            plt.savefig('%s_%s_rendered.png'%(self.series_filename.replace('.hdf5',''),corrstr),dpi=300)
            plt.show()
        

    def filter_registration(self,xshifts,yshifts,goodnesses,xmax=10,ymax=10):

        #plt.subplot(1,2,1)
        #plt.plot(xshifts)
        #plt.plot(yshifts)
        #plt.plot(goodnesses)
        #plt.subplot(1,2,2)
        #plt.plot(np.diff(xshifts))
        #plt.plot(np.diff(yshifts))
        #plt.show()
        
        xvalid = np.abs(xshifts - np.median(xshifts))<=xmax
        yvalid = np.abs(yshifts - np.median(yshifts))<=ymax
        invalid = np.logical_or(1-xvalid,1-yvalid)


        # for xv,yv,iv in zip(xvalid,yvalid,invalid):
        #     print xv,yv,iv
        # sys.exit()
        
        # plt.plot(xshifts)
        # plt.plot(yshifts)
        # plt.plot((goodnesses-np.mean(goodnesses))*10000)
        # plt.show()
        # sys.exit()

        goodnesses[invalid] = -np.inf

        return xshifts,yshifts,goodnesses
        
    def goodness_histogram(self):
        all_goodnesses = self.get_all_goodnesses()
        plt.hist(all_goodnesses,100)
        plt.show()
        
    def get_all_goodnesses(self,exclude_reference=True):
        all_goodnesses = []
        keys = self.h5.keys()
        for k in keys:
            goodnesses = self.h5['%s/%s/goodnesses'%(filename,k)][:]
            if np.prod(goodnesses)==1.0 and exclude_reference:
                continue
            all_goodnesses.append(goodnesses)
        return np.array(all_goodnesses).ravel()


    def show_images(self):

        keys = self.h5.keys()
        for k in keys:
            plt.cla()
            im = self.get_image(full_filename,vidx)
            self.imshow(im)
            plt.title(k)
            plt.pause(.1)
        plt.close()
        

if __name__=='__main__':

    frames_fn = '/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.11.21_cones/14_13_13-1T_500.hdf5' # volume 7
    h5 = H5(frames_fn)
    h5.catalog()
    
    ref = (h5['projections/ISOS'][0,:,:]+h5['projections/ISOS'][0,:,:])/2.0

    s = Series(ref,frames_fn.replace('.hdf5','')+'_registered.hdf5')

    for k in range(12):
        s.add(frames_fn,k,['ISOS','COST'],do_plot=False,overwrite=False)

    s.render(['ISOS','COST'],do_plot=True)

    
