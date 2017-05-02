import numpy as np
from matplotlib import pyplot as plt
import sys,os
from octopod import H5,utils
import glob
from scipy.ndimage import zoom
from scipy.interpolate import griddata
from scipy.signal import fftconvolve,medfilt
from scipy.io import savemat
from clicky import collector

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
        else:
            try:
                self.reference = self.h5['/reference_frame'][:,:]
            except Exception as e:
                print 'Warning: empty Series started w/o reference frame.'

    def set_reference_frame(self,reference_frame):
        self.reference = reference_frame
        self.h5.put('/reference_frame',reference_frame)

    def crop_borders(self,im):
        # remove edges of im that contain no information
        vprof = np.std(im,axis=1)
        valid_region_y = np.where(vprof)[0]
        hprof = np.std(im,axis=0)
        valid_region_x = np.where(hprof)[0]
        y1 = valid_region_y[0]
        y2 = valid_region_y[-1]
        x1 = valid_region_x[0]
        x2 = valid_region_x[-1]
        return im[y1:y2,x1:x2]

    def get_cropper(self,im):
        # return a cropper function that
        # will crop images according to im's
        # empty (uninformative) borders
        vprof = np.std(im,axis=1)
        valid_region_y = np.where(vprof)[0]
        hprof = np.std(im,axis=0)
        valid_region_x = np.where(hprof)[0]
        y1 = valid_region_y[0]
        y2 = valid_region_y[-1]
        x1 = valid_region_x[0]
        x2 = valid_region_x[-1]
        return lambda x: x[y1:y2,x1:x2]
        

    def crop_and_label_average_volume(self,tag='ISOS',labels=['ELM','ISOS','COST','RPE']):
        """Find the average volume tagged TAG and interactively
        crop it in 3D and label it"""

        av = self.h5['/average_volume/%s'%tag][:,:,:]
        zproj = av.mean(axis=1)
        xc,yc,_=collector([zproj])
        rect = np.array(xc+yc)
        rect = np.round(rect).astype(int)
        self.h5.put('average_volume/ISOS_xy_rect',rect)
        x1 = rect[0]
        x2 = rect[1]
        y1 = rect[2]
        y2 = rect[3]
        
        yproj = np.log(av.mean(axis=0)+1000.0)
        xc,yc,_=collector([yproj])
        zlim = np.round(yc).astype(int)
        self.h5.put('average_volume/ISOS_zlims',zlim)
        z1 = zlim[0]
        z2 = zlim[1]
            
        av = av[y1:y2,z1:z2,x1:x2]
        sy,sz,sx = av.shape
        
        bscan = np.log(av[sy//2,:,:]+1000.0-av.min())**.5
        for label in ['ELM','ISOS','ISOS_distal','COST_proximal','COST','RPE','OPL','CH']:
            #xc,zc,_=collector([bscan],titles=['Click upper and lower extents of %s.'%label])
            z = np.arange(bscan.shape[0])
            b = np.mean(bscan,axis=1)
            zc,xc,_=collector([(z,b)],titles=['Click left and right extents of %s.'%label])
            if len(zc)<2:
                sys.exit('Must select two points.')
            zc = np.round(np.array(zc))
            zc = np.array([zc.min(),zc.max()]).astype(int)
            self.h5.put('average_volume/%s_labels/%s'%(tag,label),zc)

    def get_labels(self,tag):
        labels = {}
        for k in self.h5['average_volume/%s_labels'%tag].keys():
            labels[k] = self.h5['average_volume/%s_labels/%s'%(tag,k)].value
        return labels
    
    def get_average_volume(self,tag):
        av = self.h5['/average_volume/%s'%tag][:,:,:]
        x1,x2,y1,y2 = self.h5['average_volume/ISOS_xy_rect'][:]
        z1,z2 = self.h5['average_volume/ISOS_zlims'][:]
        return av[y1:y2,z1:z2,x1:x2]
        
    def correct_reference_a(self,kernel_size=10,do_plot=False):
        try:
            si = self.h5['sum_image'][:,:]
            ai = self.h5['average_image'][:,:]
            ci = self.h5['counter_image'][:,:]
            eps = np.finfo(float).eps
            ci = ci + eps
            corri = self.h5['correlation_image'][:,:]

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


        corrprof = np.max(corri/ci,axis=1)[ry1:ry2]
        vprof = np.max(ci,axis=1)[ry1:ry2]
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
        if do_plot:
            clim = np.percentile(reference,(1,99.5))
            plt.figure()
            plt.imshow(reference,cmap='gray',interpolation='none',clim=clim)
            plt.title('uncorrected')
            plt.colorbar()
            plt.figure()
            plt.imshow(corrected_reference,cmap='gray',interpolation='none',clim=clim)
            plt.title('corrected a')
            plt.colorbar()
            #plt.show()


        self.h5.put('/corrected_a/reference_frame',corrected_reference)
        self.h5.put('/corrected_a/horizontal_correction',hc)
        self.h5.put('/corrected_a/vertical_correction',vc)


    def correct_reference_b(self,do_plot=False):

        files = self.h5['frames'].keys()
        frame_keys = []
        for fn in files:
                frames = self.h5['frames/%s'%fn].keys()
                for f in frames:
                    frame_keys.append('/frames/%s/%s'%(fn,f))

        all_yshifts = []
        all_xshifts = []
        for fk in frame_keys:
            ys = self.h5[fk]['y_shifts'][:]
            all_yshifts.append(ys)
            xs = self.h5[fk]['x_shifts'][:]
            all_xshifts.append(xs)

        all_yshifts = np.array(all_yshifts)
        all_xshifts = np.array(all_xshifts)

        if do_plot:
            plt.figure()
            plt.subplot(1,2,1)
            plt.imshow(all_xshifts,interpolation='none',aspect='normal')
            plt.subplot(1,2,2)
            plt.imshow(all_yshifts,interpolation='none',aspect='normal')

        # diff these:
        # we do this to effectively remove outliers,
        # since the line-to-line differences in x
        # and y can be anticipated while the absolute
        # positions cannot be (as well)
        
        d_yshifts_unfixed = np.diff(all_yshifts,axis=1)
        d_xshifts_unfixed = np.diff(all_xshifts,axis=1)

        if do_plot:
            plt.figure()
            plt.subplot(1,2,1)
            plt.imshow(d_xshifts_unfixed,interpolation='none',aspect='normal',clim=(-10,10))
            plt.title('dx, before removing outliers')
            plt.subplot(1,2,2)
            plt.imshow(d_yshifts_unfixed,interpolation='none',aspect='normal',clim=(-10,10))
            plt.title('dy, before removing outliers')

        # keep track of the medians of the first columns,
        # because we need these values to reconstruct
        # the absolute positions with cumsum below
        first_y = np.median(all_yshifts[:,0])
        first_x = np.median(all_xshifts[:,0])
        # remove outliers (replace with nans)
        def outlier_to_nan(arr,ulim=.5,llim=-.5):
            invalid = np.where(np.logical_or(arr>ulim,arr<llim))
            arr[invalid] = np.nan
            return arr

        d_yshifts = outlier_to_nan(d_yshifts_unfixed)
        d_xshifts = outlier_to_nan(d_xshifts_unfixed)

        if do_plot:
            plt.figure()
            plt.subplot(1,2,1)
            plt.title('dx, after removing outliers')
            plt.imshow(d_xshifts,interpolation='none',aspect='normal',clim=(-10,10))
            plt.subplot(1,2,2)
            plt.imshow(d_yshifts,interpolation='none',aspect='normal',clim=(-10,10))
            plt.title('dy, after removing outliers')
            
        # now, nanmean these diffs to get the average
        # lag biases; these will be integrated to get
        # the eye position
        d_y = np.nanmean(d_yshifts,axis=0)
        d_x = np.nanmean(d_xshifts,axis=0)
        d_y = np.array([first_y]+list(d_y))
        d_x = np.array([first_x]+list(d_x))

        vc = np.cumsum(d_y)
        hc = np.cumsum(d_x)

        if False: # save eye motion data to a MAT file
            outdict = {}
            outdict['horizontal_eye_position_pixels'] = hc
            outdict['vertical_eye_position_pixels'] = vc
            savemat('eye_motion.mat',outdict)

        reference = self.h5['reference_frame'][:,:]
        original_sy,original_sx = reference.shape
        
        # downsample the movement estimates:
        for fk in self.h5['frames'].keys():
            for idx in self.h5['frames'][fk].keys():
                oversample_factor = self.h5['frames'][fk][idx]['oversample_factor'].value
                break

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

        clim = np.percentile(reference,(1,99.5))
        
        if do_plot:
            plt.figure()
            plt.imshow(reference,cmap='gray',interpolation='none',clim=clim)
            plt.title('uncorrected')
            plt.colorbar()
            plt.figure()
            plt.imshow(corrected_reference,cmap='gray',interpolation='none',clim=clim)
            plt.title('corrected b')
            plt.colorbar()
            plt.show()

        self.h5.put('/corrected_b/reference_frame',corrected_reference)
        self.h5.put('/corrected_b/horizontal_correction',hc)
        self.h5.put('/corrected_b/vertical_correction',vc)
        
        
        
    def add(self,filename,vidx,layer_names=None,overwrite=False,oversample_factor=3,strip_width=3.0,do_plot=False,use_gaussian=True):
        
        print 'Adding %s, volume %d.'%(filename,vidx)
        
        target_tag = self.tag_template%(os.path.split(filename)[1],vidx)

        if self.h5.has('/frames/%s'%target_tag) and not overwrite:
            print 'Series already has entry for %s.'%target_tag
            return

        target,label = self.get_image(filename,vidx,layer_names)
        reference = self.reference
        y,x,g = utils.strip_register(target,reference,oversample_factor,strip_width,do_plot=do_plot,use_gaussian=use_gaussian)
        
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

        label = '_'.join(layer_names)

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
        return out,label

    def get_volume(self,filename_stub,vidx,data_block):
        filename = os.path.join(self.working_directory,filename_stub)
        target_h5 = H5(filename)
        return target_h5[data_block][vidx,:,:,:]

    def get_z_offsets(self,filename_stub,vidx,layer_name=None):
        filename = os.path.join(self.working_directory,filename_stub)
        target_h5 = H5(filename)
        if layer_name is None:
            out = target_h5['/model/z_offsets'][vidx,:,:]
        else:
            out = target_h5['/model/volume_labels/%s'%layer_name][vidx,:,:]
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

            
    def is_volume_rendered(self):
        out = True
        try:
            test = self.h5['/sum_volume']
        except:
            out = False
        return out

            
    def is_corrected_a(self):
        out = True
        try:
            test = self.h5['/corrected_a']
        except:
            out = False

        return out

    
    def is_corrected_b(self):
        out = True
        try:
            test = self.h5['/corrected_b']
        except:
            out = False

        return out

    def goodness_histogram(self):
        files = self.h5['frames'].keys()
        all_goodnesses = []
        for filename in files:
            keys = self.h5['/frames/%s'%filename].keys()
            for k in keys:
                print filename
                goodnesses = list(self.h5['/frames/%s/%s/goodnesses'%(filename,k)][:])
                all_goodnesses = all_goodnesses + goodnesses

        plt.hist(all_goodnesses,500)
        plt.show()
    
    def render(self,layer_names=None,goodness_threshold=0.0,correlation_threshold=-1.0,overwrite=False,oversample_factor=3,do_plot=False,frames_to_save=[],left_crop=0):

        if len(frames_to_save):
            frames_directory = self.series_filename.replace('.hdf5','')+'_saved_frames'
            if not os.path.exists(frames_directory):
                os.makedirs(frames_directory)

        files = self.h5['frames'].keys()
        
        sign = -1
        # remember the convention here: x and y shifts are the
        # amount of shift required to align the line in question
        # with the reference image
        # first, find the minimum and maximum x and y shifts,
        # in order to know how big the canvas must be for fitting
        # all of the lines
        xmin = np.inf
        xmax = -np.inf
        ymin = np.inf
        ymax = -np.inf

        reg_dict = {}
        
        for filename in files:
            keys = self.h5['/frames/%s'%filename].keys()
            for k in keys:
                test,label = self.get_image(filename,0,layer_names)
                n_slow,n_fast = test.shape

                goodnesses = self.h5['/frames/%s/%s/goodnesses'%(filename,k)][:]
                xshifts = sign*self.h5['/frames/%s/%s/x_shifts'%(filename,k)][:]
                yshifts = sign*self.h5['/frames/%s/%s/y_shifts'%(filename,k)][:]

                xshifts = np.squeeze(xshifts)
                yshifts = np.squeeze(yshifts)

                xshifts,yshifts,goodnesses,valid = self.filter_registration(xshifts,yshifts,goodnesses)

                newxmin = np.min(xshifts)
                newxmax = np.max(xshifts)
                newymin = np.min(yshifts)
                newymax = np.max(yshifts)
                
                xmin = min(xmin,newxmin)
                xmax = max(xmax,newxmax)
                ymin = min(ymin,newymin)
                ymax = max(ymax,newymax)

                yshifts = yshifts + valid

                reg_dict[(filename,k)] = (xshifts,yshifts,goodnesses,valid)


        canvas_width = xmax-xmin+n_fast
        canvas_height = ymax-ymin+n_slow

        ref_x1 = 0 - xmin
        ref_y1 = 0 - ymin
        ref_x2 = ref_x1 + n_fast
        ref_y2 = ref_y1 + n_slow
        
        self.h5.put('/reference_coordinates/x1',ref_x1)
        self.h5.put('/reference_coordinates/x2',ref_x2)
        self.h5.put('/reference_coordinates/y1',ref_y1)
        self.h5.put('/reference_coordinates/y2',ref_y2)
        
        for key in reg_dict.keys():
            xs,ys,g,v = reg_dict[key]
            xs = xs - xmin
            ys = ys - ymin
            reg_dict[key] = (xs,ys,g,v)
                

        canvas_width = canvas_width*oversample_factor
        canvas_height = (canvas_height+1)*oversample_factor

        rmean = np.mean(self.reference)
        embedded_reference = np.ones((canvas_height,canvas_width))*rmean
        ref_x1 = ref_x1*oversample_factor
        ref_x2 = ref_x2*oversample_factor
        ref_y1 = ref_y1*oversample_factor
        ref_y2 = ref_y2*oversample_factor
        ref_oversampled = zoom(self.reference,oversample_factor)
        embedded_reference[ref_y1:ref_y2,ref_x1:ref_x2] = ref_oversampled
        
        sum_image = np.zeros((canvas_height,canvas_width))
        counter_image = np.zeros((canvas_height,canvas_width))
        correlation_image = np.zeros((canvas_height,canvas_width))

        for k in reg_dict.keys():
            xshifts,yshifts,goodnesses,indices = reg_dict[k]
            filename = k[0]
            frame_index = int(k[1])
            im,label = self.get_image(filename,frame_index,layer_names)
            correlation_vector = np.zeros((n_slow))
            
            for idx,xs,ys,g in zip(indices,xshifts,yshifts,goodnesses):
                xs = xs + left_crop
                line = im[idx,left_crop:]
                line = np.expand_dims(line,0)
                block = zoom(line,oversample_factor)
                bsy,bsx = block.shape
                x1 = int(np.round(xs*oversample_factor))
                x2 = x1 + bsx
                y1 = int(np.round(ys*oversample_factor))
                y2 = y1 + bsy

                ref_section = embedded_reference[y1:y2,x1:x2].ravel()
                #eref_std = np.std(ref_section)
                corr = np.corrcoef(ref_section,block.ravel())[1,0]
                correlation_image[y1:y2,x1:x2] = correlation_image[y1:y2,x1:x2] + corr
                correlation_vector[idx] = corr
                if corr>correlation_threshold:
                    sum_image[y1:y2,x1:x2] = sum_image[y1:y2,x1:x2] + block
                    counter_image[y1:y2,x1:x2] = counter_image[y1:y2,x1:x2] + 1.0
            self.h5.put('/frames/%s/%s/correlations'%(filename,k[1]),correlation_vector)
            print correlation_vector.max(),correlation_vector.min(),correlation_vector.mean()
            if do_plot:
                temp = counter_image.copy()
                temp[np.where(temp==0)] = 1.0
                av = sum_image/temp
                plt.clf()
                plt.subplot(1,3,1)
                plt.cla()
                plt.imshow(av,cmap='gray',interpolation='none')

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

            
        temp = counter_image.copy()
        temp[np.where(temp==0)] = 1.0
        av = sum_image/temp

        cropper = self.get_cropper(counter_image)
        self.h5.put('/correlation_image/%s'%label,cropper(correlation_image))
        self.h5.put('/counter_image/%s'%label,cropper(counter_image))
        self.h5.put('/average_image/%s'%label,cropper(av))
        self.h5.put('/sum_image/%s'%label,cropper(sum_image))
        
        if do_plot:
            plt.close()

            plt.subplot(1,2,1)
            plt.imshow(av,cmap='gray',interpolation='none')

            plt.subplot(1,2,2)
            plt.imshow(counter_image)
            plt.colorbar()

            #plt.savefig('%s_%s_rendered.png'%(self.series_filename.replace('.hdf5',''),corrstr),dpi=300)
            plt.show()
        
            
    def render_volume(self,layer_names=None,goodness_threshold=0.0,correlation_threshold=-1.0,overwrite=False,oversample_factor=3,do_plot=False,left_crop=0,data_block='/flattened_data',offset_medfilt_kernel=9,layer_name='ISOS',align_bscan=False):

        files = self.h5['frames'].keys()
        print files
        
        sign = -1
        # remember the convention here: x and y shifts are the
        # amount of shift required to align the line in question
        # with the reference image
        # first, find the minimum and maximum x and y shifts,
        # in order to know how big the canvas must be for fitting
        # all of the lines
        xmin = np.inf
        xmax = -np.inf
        ymin = np.inf
        ymax = -np.inf
        zmin = np.inf
        zmax = -np.inf
        max_depth = -np.inf

        reg_dict = {}

        def fix_corners(im):
            fill_value = np.median(im)
            to_check = np.zeros(im.shape)
            rad = (offset_medfilt_kernel-1)//2
            to_check[:rad,:rad] = 1
            to_check[-rad:,:rad] = 1
            to_check[:rad,-rad:] = 1
            to_check[-rad:,-rad:] = 1
            im[np.where(np.logical_and(to_check,1-im))] = fill_value
            return im
        
        for filename in files:
            keys = self.h5['/frames/%s'%filename].keys()
            for k in keys:
                test = self.get_volume(filename,int(k),data_block)
                test = np.abs(test[:,:,left_crop:])
                orig_vol_shape = test.shape
                zshifts = self.get_z_offsets(filename,int(k),layer_name)[:,left_crop:]
                zshifts = medfilt(zshifts,offset_medfilt_kernel)
                zshifts = fix_corners(zshifts)
                
                n_slow,n_depth,n_fast = test.shape
                if n_depth>max_depth:
                    max_depth = n_depth
                if False:
                    for k in range(n_depth):
                        im = test[:,k,:]
                        plt.cla()
                        plt.imshow(im,cmap='gray',interpolation='none',clim=np.percentile(test,(5,99.95)))
                        plt.pause(.1)
                    sys.exit()
                
                goodnesses = self.h5['/frames/%s/%s/goodnesses'%(filename,k)][:]
                xshifts = sign*self.h5['/frames/%s/%s/x_shifts'%(filename,k)][:]
                yshifts = sign*self.h5['/frames/%s/%s/y_shifts'%(filename,k)][:]

                xshifts = np.squeeze(xshifts)
                yshifts = np.squeeze(yshifts)

                xshifts,yshifts,goodnesses,valid = self.filter_registration(xshifts,yshifts,goodnesses)

                zshifts = zshifts[valid,:]
                
                newxmin = np.min(xshifts)
                newxmax = np.max(xshifts)
                newymin = np.min(yshifts)
                newymax = np.max(yshifts)
                newzmin = np.min(zshifts)
                newzmax = np.max(zshifts)

                
                xmin = min(xmin,newxmin)
                xmax = max(xmax,newxmax)
                ymin = min(ymin,newymin)
                ymax = max(ymax,newymax)
                zmin = min(zmin,newzmin)
                zmax = max(zmax,newzmax)
                yshifts = yshifts + valid

                reg_dict[(filename,k)] = (xshifts,yshifts,zshifts,goodnesses,valid)

        
        canvas_width = xmax-xmin+n_fast
        canvas_height = ymax-ymin+n_slow
        canvas_depth = zmax-zmin+max_depth+1
        print 'canvas_depth start',canvas_depth
        
        ref_x1 = 0 - xmin
        ref_y1 = 0 - ymin
        ref_x2 = ref_x1 + n_fast
        ref_y2 = ref_y1 + n_slow
        
        self.h5.put('/reference_coordinates/x1',ref_x1)
        self.h5.put('/reference_coordinates/x2',ref_x2)
        self.h5.put('/reference_coordinates/y1',ref_y1)
        self.h5.put('/reference_coordinates/y2',ref_y2)
        
        for key in reg_dict.keys():
            xs,ys,zs,g,v = reg_dict[key]
            xs = xs - xmin
            ys = ys - ymin
            zs = zs - zmin
            reg_dict[key] = (xs,ys,zs,g,v)

        if False:
            for key in reg_dict.keys():
                plt.clf()
                plt.imshow(reg_dict[key][2])
                plt.colorbar()
                plt.pause(1)
            sys.exit()
            
        canvas_width = canvas_width*oversample_factor
        canvas_height = (canvas_height+1)*oversample_factor
        canvas_depth0 = canvas_depth
        canvas_depth = canvas_depth*oversample_factor
        print 'canvas_depth oversampled',canvas_depth
        sum_image = np.zeros((canvas_height,canvas_depth,canvas_width))
        counter_image = np.ones((canvas_height,canvas_depth,canvas_width))*1e-10

        
        for volume_count,k in enumerate(reg_dict.keys()):
            xshifts,yshifts,zshifts,goodnesses,indices = reg_dict[k]
            filename = k[0]
            frame_index = int(k[1])
            vol = np.abs(self.get_volume(filename,frame_index,data_block))[:,:,left_crop:]
            
            for bscan_count,(idx,xs,ys,zs,g) in enumerate(zip(indices,xshifts,yshifts,zshifts,goodnesses)):
                bscan = vol[idx,:,:]
                n_depth,n_fast = bscan.shape
                shifted_bscan = np.zeros((canvas_depth0,n_fast))
                print 'bscan shape',bscan.shape
                print 'shifted_bscan shape',shifted_bscan.shape
                print 'zs lims',zs.min(),zs.max()
                print 'zmax',zmax
                if align_bscan:
                    for i_fast in range(n_fast):
                        z1 = zmax-zs[i_fast]
                        z2 = z1 + n_depth
                        shifted_bscan[z1:z2,i_fast] = bscan[:,i_fast]
                else:
                    z1 = zmax - int(round(np.mean(zs)))
                    z2 = z1 + n_depth
                    cut_count = 0
                    while z2>shifted_bscan.shape[0]:
                        bscan = bscan[:-1,:]
                        z2 = z2 - 1
                        cut_count =+ 1
                    if cut_count:
                        print 'cut %d lines'%cut_count    
                    shifted_bscan[z1:z2,:] = bscan
                print

                if False:
                    plt.subplot(3,1,1)
                    plt.imshow(bscan,interpolation='none',cmap='gray',aspect='normal')
                    plt.subplot(3,1,2)
                    plt.imshow(shifted_bscan,interpolation='none',cmap='gray',aspect='normal')
                    plt.subplot(3,1,3)
                    plt.plot(zs)
                    plt.figure()
                    plt.plot(np.mean(bscan,axis=1))
                    plt.plot(np.mean(shifted_bscan,axis=1))
                    plt.show()
                    continue
                
                block = np.zeros((oversample_factor,canvas_depth,n_fast*oversample_factor))
                shifted_bscan = zoom(shifted_bscan,oversample_factor)
                
                for ofk in range(oversample_factor):
                    block[ofk,:,:] = shifted_bscan
                bsy,bsz,bsx = block.shape
                
                x1 = int(np.round(xs*oversample_factor))
                x2 = x1 + bsx
                y1 = int(np.round(ys*oversample_factor))
                y2 = y1 + bsy
                
                sum_image[y1:y2,:,x1:x2] = sum_image[y1:y2,:,x1:x2] + block
                counter_image[y1:y2,:,x1:x2] = counter_image[y1:y2,:,x1:x2] + 1.0
                print volume_count,bscan_count
            if do_plot:
                vav = sum_image[:,500:600,200]/counter_image[:,500:600,200]
                hav = sum_image[200,500:600,:]/counter_image[200,500:600,:]
                plt.clf()
                plt.subplot(1,2,1)
                plt.cla()
                plt.imshow(vav.T,cmap='gray',interpolation='none',aspect='normal')
                plt.colorbar()
                plt.title(volume_count)
                plt.subplot(1,2,2)
                plt.cla()
                plt.imshow(hav,cmap='gray',interpolation='none',aspect='normal')
                plt.colorbar()
                plt.pause(.0000000001)

            
        av = sum_image/counter_image
        label = layer_name
        self.h5.put('/counter_volume/%s'%label,counter_image)
        self.h5.put('/average_volume/%s'%label,av)
        self.h5.put('/sum_volume/%s'%label,sum_image)
        
        if do_plot:
            plt.close()
            plt.subplot(1,2,1)
            plt.cla()
            plt.imshow(vav.T,cmap='gray',interpolation='none',aspect='normal')
            plt.colorbar()
            plt.title(volume_count)
            plt.subplot(1,2,2)
            plt.cla()
            plt.imshow(hav,cmap='gray',interpolation='none',aspect='normal')
            plt.colorbar()
            plt.show()
        
            
    def filter_registration(self,xshifts,yshifts,goodnesses,xmax=3,ymax=3,medfilt_region=49,do_plot=False):
        xmed = medfilt(xshifts,medfilt_region)
        ymed = medfilt(yshifts,medfilt_region)
        xerr = np.abs(xshifts-xmed)
        yerr = np.abs(yshifts-ymed)
        xvalid = xerr<=xmax
        yvalid = yerr<=ymax
        

        valid = np.where(np.logical_and(xvalid,yvalid))[0]

        #print '%d points: '%len(valid),valid
        if do_plot:
            plt.figure()
            plt.subplot(1,2,1)
            plt.plot(xshifts,'k-')
            #plt.plot(xmed,'b-')
            plt.plot(yshifts,'k--')
            #plt.plot(ymed,'b--')

            plt.subplot(1,2,2)
            plt.plot(valid,xshifts[valid],'ks')
            plt.plot(valid,yshifts[valid],'ko')
            
            plt.show()

        return xshifts[valid],yshifts[valid],goodnesses[valid],valid
        
        


    def show_images(self):

        keys = self.h5.keys()
        for k in keys:
            plt.cla()
            im,label = self.get_image(full_filename,vidx)
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

    
