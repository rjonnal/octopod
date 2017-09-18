import numpy as np
from matplotlib import pyplot as plt
import sys,os
from octopod import H5,utils
import glob
from scipy.ndimage import zoom
from scipy.interpolate import griddata
from scipy.signal import fftconvolve,medfilt
from scipy.optimize import curve_fit
from scipy.io import savemat
from clicky import collector
import matplotlib.colors as colors
import matplotlib.cm as cmx



class Cone:

    def __init__(self,vol,isos_idx,cost_idx,x,y,origin='',x0=None,y0=None,properties={}):
        """This class stores information about a single cone. Required
        parameters are a 3D array representing the cone's reflectivity;
        indices of the ISOS and COST reflections; the x and y
        coordinates of the cone in the reference coordinate
        space. Optional params are a string describing the cone's
        original data set (e.g. the dataset tag and volume index); x0
        and y0, the coordinates of the cone in its original volume; and
        a dictionary of other properties, e.g. whether the cone had been
        stimulated or not, its spectral class, etc."""
        
        self.iidx = isos_idx
        self.cidx = cost_idx
        self.vol = vol
        self.x = x
        self.y = y
        self.origin = origin
        self.x0 = x0
        self.y0 = y0
        self.properties = properties
        self.prof = np.mean(np.mean(np.abs(vol),1),0)

    def gaussian_mixture_model(self,do_plot=False):

        if do_plot:
            plt.figure()
        def gmm(x,x0,x1,s0,s1,a0,a1):
            XX0 = np.arange(len(x)) - x0
            XX1 = np.arange(len(x)) - x1
            y0 = np.exp(-(XX0**2)/(2*s0**2))*a0
            y1 = np.exp(-(XX1**2)/(2*s1**2))*a1
            if do_plot:
                plt.cla()
                plt.plot(x,y0+y1)
                plt.plot(x,self.prof)
                plt.pause(.0001)
            return y0+y1

        # initial guess:
        p_a0 = self.prof[self.iidx]
        p_a1 = self.prof[self.cidx]
        p_s0 = 1.0
        p_s1 = 1.0
        p_x0 = self.iidx
        p_x1 = self.cidx
        p = [p_x0,p_x1,p_s0,p_s1,p_a0,p_a1]
        lower = [p_x0-1.0,p_x1-1.0,0.1,0.1,p_a0*.75,p_a1*.75]
        upper = [p_x0+1.0,p_x1+1.0,3.0,3.0,p_a0*1.25,p_a1*1.25]
        #bounds = (lower,upper)
        
        aa,bb = curve_fit(gmm,np.arange(len(self.prof)),self.prof,p0=p)
        
        if do_plot:
            plt.close()

        print 'bye'
        return aa,bb
        
        


class Series:

    def __init__(self,series_filename,reference_frame=None):
        self.h5 = H5(series_filename)
        self.tag = os.path.splitext(series_filename)[0]
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

    def make_cone_catalog(self,points,minimum_goodness=10.0,output_radius=2,match_radius=2.0,do_plot=False):
        """Take a set of points in this Series' reference image,
        corresponding to the x and y coordinate of cone centers,
        and identify and crop the corresponding cone out of all
        this Series' volumes."""
        
        fkeys = self.h5['/frames'].keys()
        for fkidx,fk in enumerate(fkeys):
            dataset_fn = os.path.join(self.working_directory,fk)
            dataset_h5 = H5(dataset_fn)
            ikeys = self.h5['/frames'][fk].keys()
            for ikidx,ik in enumerate(ikeys):
                vol = dataset_h5['/flattened_data'][int(ik),:,:,:]
                cost_depths = dataset_h5['model/volume_labels/COST'][int(ik),:,:]
                isos_depths = dataset_h5['model/volume_labels/ISOS'][int(ik),:,:]

                x = self.h5['/frames'][fk][ik]['x_shifts'][:]
                y = self.h5['/frames'][fk][ik]['y_shifts'][:]
                c = self.h5['/frames'][fk][ik]['correlations'][:]
                g = self.h5['/frames'][fk][ik]['goodnesses'][:]

                model_profile = dataset_h5['model/profile'][:]
                model_isos = dataset_h5['model/labels/ISOS'].value
                model_cost = dataset_h5['model/labels/COST'].value
                yramp = np.arange(len(y)) - y
                # yramp[n] is now the location of the nth target row in the reference space
                # so we need to find n such that yramp[n]-pty is minimized
                
                cmed = int(np.median(cost_depths))
                imed = int(np.median(isos_depths))
                volume_enface_projection = np.abs(vol[:,imed-2:cmed+2,:]).mean(axis=1)

                border = 3
                
                for idx,(ptx,pty) in enumerate(points):
                    print 'frame %d/%d; vol %d/%d; cone %d/%d at %d,%d'%(fkidx+1,len(fkeys),ikidx+1,len(ikeys),idx+1,len(points),ptx,pty)
                    yerr = np.abs(pty-yramp)
                    match_index = np.argmin(yerr)
                    if yerr[match_index]<=match_radius and g[match_index]>minimum_goodness:
                        print 'match exists'
                        # get the target coordinates, and then ascend the target en face projection
                        # to the peak (center) of the cone:
                        yout = int(match_index)
                        xout = int(ptx + x[match_index])
                        xout,yout = utils.ascend2d(volume_enface_projection,xout,yout,do_plot=False)
                        
                        # 3D-crop the cone out of the volume, and make an axial profile from it
                        try:
                            cone_volume = self.get_subvol(vol,yout,0,xout,output_radius,np.inf,output_radius)
                            cone_profile = np.mean(np.mean(np.abs(cone_volume),axis=2),axis=0)
                            shift,corr = utils.nxcorr(model_profile,cone_profile)
                        except Exception as e:
                            print e
                            continue

                        # get some overall stats of the volume
                        # we'll write these to the hdf5 file later
                        avol = np.abs(cone_volume)
                        aprof = np.mean(np.mean(avol,axis=2),axis=0)
                        noise_floor = aprof[np.argsort(aprof)[:10]]
                        noise_mean = noise_floor.mean()
                        noise_std = noise_floor.std()
                        full_volume_mean = avol.mean()
                        full_volume_std = avol.std()
                        full_profile_mean = aprof.mean()
                        full_profile_std = aprof.std()
                        
                        
                        # let's try to label this cone's peaks
                        #print 'mp',model_profile
                        #print 'cp',cone_profile
                        isos_guess = int(model_isos + shift)
                        cost_guess = int(model_cost + shift)

                        # now use the height of the cone_profile at
                        # these guesses, combined with locations and
                        # distances of other peaks to refine the guesses
                        max_displacement = 4
                        peaks = utils.find_peaks(cone_profile)
                        heights = cone_profile[peaks]/cone_profile.std()
                        isos_dpeaks = np.abs(peaks-isos_guess)
                        isos_dpeaks[np.where(isos_dpeaks>=max_displacement)] = 2**16
                        cost_dpeaks = np.abs(peaks-cost_guess)
                        cost_dpeaks[np.where(cost_dpeaks>=max_displacement)] = 2**16
                        isos_scores = heights - isos_dpeaks
                        cost_scores = heights - cost_dpeaks
                        isos_guess = peaks[np.argmax(isos_scores)]
                        cost_guess = peaks[np.argmax(cost_scores)]

                        # now cross correlate the cuts through the cone with
                        # cone_profile to see if axial eye movement between
                        # b-scans has resulted in a peak shift
                        sy,sz,sx = cone_volume.shape
                        colors = 'rgbkcym'
                        sheet = []
                        for idx,coney in enumerate(range(sy)):
                            color = colors[idx%len(colors)]
                            cut_profile = np.abs(cone_volume[coney,:,:]).mean(axis=1)
                            shift,corr = utils.nxcorr(cut_profile,cone_profile)
                            z1 = int(isos_guess-border-shift)
                            z2 = int(cost_guess+border-shift)
                            cut = cone_volume[coney,z1:z2,:]
                            sheet.append(cut)

                        sheet = np.array(sheet)
                        sheet = np.transpose(sheet,(0,2,1))
                        point_string = '%d_%d'%(ptx,pty)
                        os_length = cost_guess-isos_guess
                        key_root = '/cone_catalog/%s/%s/%s'%(point_string,fk,ik)


                        c = Cone(sheet,border,border+os_length,ptx,pty)

                        try:
                            fit,fit_covar = c.gaussian_mixture_model(True)
                        except Exception as e:
                            print e
                            fit,fit_covar = c.gaussian_mixture_model(True)
                            
                        fit_isos_z,fit_cost_z,fit_isos_sigma,fit_cost_sigma,fit_isos_amplitude,fit_cost_amplitude = fit
                        
                        perr = np.sqrt(np.diag(fit_covar))
                        fit_isos_z_error,fit_cost_z_error,fit_isos_sigma_error,fit_cost_sigma_error,fit_isos_amplitude_error,fit_cost_amplitude_error = perr

                        self.h5.put('%s/x'%key_root,xout)
                        print 'hi'
                        self.h5.put('%s/y'%key_root,yout)
                        print 'hi'
                        self.h5.put('%s/isos_z'%key_root,border)
                        print 'hi'
                        self.h5.put('%s/cost_z'%key_root,border+os_length)
                        self.h5.put('%s/cone_volume'%key_root,sheet)
                        self.h5.put('%s/noise_mean'%key_root,noise_mean)
                        self.h5.put('%s/noise_std'%key_root,noise_std)
                        self.h5.put('%s/full_volume_mean'%key_root,full_volume_mean)
                        self.h5.put('%s/full_volume_std'%key_root,full_volume_std)
                        self.h5.put('%s/full_profile_mean'%key_root,full_profile_mean)
                        self.h5.put('%s/full_profile_std'%key_root,full_profile_std)
                        self.h5.put('%s/fit/isos_z'%key_root,fit_isos_z)
                        self.h5.put('%s/fit/cost_z'%key_root,fit_cost_z)
                        self.h5.put('%s/fit/isos_sigma'%key_root,fit_isos_sigma)
                        self.h5.put('%s/fit/cost_sigma'%key_root,fit_cost_sigma)
                        self.h5.put('%s/fit/isos_amplitude'%key_root,fit_isos_amplitude)
                        self.h5.put('%s/fit/cost_amplitude'%key_root,fit_cost_amplitude)
                        self.h5.put('%s/fit/isos_z_error'%key_root,fit_isos_z_error)
                        self.h5.put('%s/fit/cost_z_error'%key_root,fit_cost_z_error)
                        self.h5.put('%s/fit/isos_sigma_error'%key_root,fit_isos_sigma_error)
                        self.h5.put('%s/fit/cost_sigma_error'%key_root,fit_cost_sigma_error)
                        self.h5.put('%s/fit/isos_amplitude_error'%key_root,fit_isos_amplitude_error)
                        self.h5.put('%s/fit/cost_amplitude_error'%key_root,fit_cost_amplitude_error)

                        
    def get_n_volumes(self):
        count = 0
        fkeys = self.h5['/frames'].keys()
        for fk in fkeys:
            ikeys = self.h5['/frames'][fk].keys()
            count = count + len(ikeys)
        return count


    def get_volume_dictionary(self,order=None):
        d = {}
        counter = 0
        frame_keys = self.h5['/frames'].keys()
        if order is not None:
            frame_keys = frame_keys[order]
        for fk in frame_keys:
            index_keys = self.h5['/frames/%s'%fk].keys()
            for ik in index_keys:
                d[(fk,ik)] = counter
                counter = counter + 1
        return d
                
    
    def get_n_cones(self):
        try:
            n_cones = self.h5['cone_catalog/globals/n_cones'].value
        except Exception as e:
            cone_catalog = self.h5['cone_catalog']
            cone_keys = cone_catalog.keys()
            n_cones = len(cone_keys)
            self.h5.put('/cone_catalog/globals/n_cones',n_cones)
        return n_cones
    
    def get_cone_volume_size(self):
        fast_maxes = []
        try:
            slow_max = self.h5['cone_catalog/globals/slow_max'].value
            fast_max = self.h5['cone_catalog/globals/fast_max'].value
            depth_max = self.h5['cone_catalog/globals/depth_max'].value
        except Exception as e:
            cone_catalog = self.h5['cone_catalog']
            cone_keys = cone_catalog.keys()
            slow_max,fast_max,depth_max = 0,0,0
            for ck in cone_keys:
                if ck in ['globals']:
                    continue
                frame_keys = cone_catalog['%s'%ck].keys()
                for fk in frame_keys:
                    index_keys = cone_catalog['%s/%s'%(ck,fk)].keys()
                    for ik in index_keys:
                        dims = cone_catalog['%s/%s/%s/cone_volume'%(ck,fk,ik)].shape
                        print ck,fk,dims
                        if dims[1]>5:
                            continue
                        if dims[0]>slow_max:
                            slow_max = dims[0]
                        if dims[1]>fast_max:
                            fast_max = dims[1]
                        if dims[2]>depth_max:
                            depth_max = dims[2]
                        fast_maxes.append(dims[1])
                        
            self.h5.put('/cone_catalog/globals/slow_max',slow_max)
            self.h5.put('/cone_catalog/globals/fast_max',fast_max)
            self.h5.put('/cone_catalog/globals/depth_max',depth_max)
        return slow_max,fast_max,depth_max
                    
        

    def analyze_cone_phase(self,do_plot=True):
        av = self.crop_average_to_reference()

        clim = np.percentile(av,(1,99.5))
        
        if do_plot:
            jet = cm = plt.get_cmap('jet') 
            cNorm  = colors.Normalize(vmin=-1.0, vmax=1.0)
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
            sy,sx = av.shape
            ar = float(sx)/float(sy)
            dpi = 100.0
            outx = 800.0/dpi
            outy = outx/ar
            plt.figure(2)
            plt.figure(1,figsize=(outx,outy))
            plt.axes([0,0,1,1])
            plt.imshow(av,clim=clim,cmap='gray',interpolation='none')
            plt.autoscale(False)
            outfn = '%s_unmarked_average.png'%self.tag
            plt.savefig(outfn)
            
        volume_dictionary = self.get_volume_dictionary()
        
        n_volumes = self.get_n_volumes()
        n_cones = self.get_n_cones()

        cone_catalog = self.h5['cone_catalog']
        cone_keys = cone_catalog.keys()

        
        for cone_index,ck in enumerate(cone_keys):
            stim_d_phase = []
            no_stim_d_phase = []
            if ck in ['slow_max','fast_max','depth_max','n_cones']:
                continue
            frame_keys = cone_catalog['%s'%ck].keys()

            coords = ck.split('_')
            x0 = int(coords[0])
            y0 = int(coords[1])
            
            for fk in frame_keys:
                if fk.find('no_stimulus')>-1:
                    stimulus = 0
                elif fk.find('stimulus')>-1:
                    stimulus = 1
                else:
                    stimulus = 0
                    
                index_keys = cone_catalog['%s/%s'%(ck,fk)].keys()
                for ik in index_keys:
                    vol = cone_catalog['%s/%s/%s/cone_volume'%(ck,fk,ik)][:,:,:]
                    
                    isos_depth = cone_catalog['%s/%s/%s/isos_z'%(ck,fk,ik)].value
                    cost_depth = cone_catalog['%s/%s/%s/cost_z'%(ck,fk,ik)].value

                    x = cone_catalog['%s/%s/%s/x'%(ck,fk,ik)].value
                    y = cone_catalog['%s/%s/%s/y'%(ck,fk,ik)].value
                    cone_origin = '%s_%s'%(fk,ik)
                    xy0 = ck.split('_')
                    x0 = int(xy0[0])
                    y0 = int(xy0[1])

                    c = Cone(vol,isos_depth,cost_depth,x,y,cone_origin,x0,y0)
                    
                    isos = vol[:,:,isos_depth].ravel()
                    cost = vol[:,:,cost_depth].ravel()
                    
                    
                    isos_phase = np.angle(isos)
                    cost_phase = np.angle(cost)
                    isos_amp = np.abs(isos)
                    cost_amp = np.abs(cost)
                    isos_real = np.real(isos)
                    cost_real = np.real(cost)
                    isos_imag = np.imag(isos)
                    cost_imag = np.imag(cost)

                    
                    amp = (isos_amp+cost_amp)/2.0

                    valid = np.where(amp>amp.mean()-0*amp.std())[0]
                    d_phase = (cost_phase-isos_phase)
                    d_phase = utils.unwrap(d_phase)
                    d_phase_median = np.median(d_phase[valid])

                    
                    if stimulus:
                        stim_d_phase.append(d_phase_median)
                    else:
                        no_stim_d_phase.append(d_phase_median)

                    
                    #print x0,y0,x,y,fk,ik,stimulus,d_phase_median
            pv = np.nanvar(stim_d_phase)/np.nanvar(no_stim_d_phase)
            print x0,y0,np.var(no_stim_d_phase),np.nanvar(stim_d_phase),pv
            if do_plot and pv>0.01:
                plt.figure(1)
                lu = min(1.0,-np.log(np.var(stim_d_phase)))
                colorVal = scalarMap.to_rgba(lu)
                plt.plot(x0*3,y0*3,marker='o',color=colorVal,alpha=0.5,markersize=10)
                
        if do_plot:
            outfn = '%s_marked_average.png'%self.tag
            plt.figure(1)
            plt.savefig(outfn)
            plt.show()
            
    def fit_cones(self):
        
        cone_catalog = self.h5['cone_catalog']
        cone_keys = cone_catalog.keys()
        
        for cone_index,ck in enumerate(cone_keys):
            if ck in ['slow_max','fast_max','depth_max','n_cones']:
                continue
            frame_keys = cone_catalog['%s'%ck].keys()

            coords = ck.split('_')
            x = int(coords[0])
            y = int(coords[1])
            
            for fk in frame_keys:
                    
                index_keys = cone_catalog['%s/%s'%(ck,fk)].keys()
                for ik in index_keys:
                    vol = cone_catalog['%s/%s/%s/cone_volume'%(ck,fk,ik)][:,:,:]
                    
                    isos_depth = cone_catalog['%s/%s/%s/isos_z'%(ck,fk,ik)].value
                    cost_depth = cone_catalog['%s/%s/%s/cost_z'%(ck,fk,ik)].value

                    x0 = cone_catalog['%s/%s/%s/x'%(ck,fk,ik)].value
                    y0 = cone_catalog['%s/%s/%s/y'%(ck,fk,ik)].value
                    cone_origin = '%s_%s'%(fk,ik)

                    c = Cone(vol,isos_depth,cost_depth,x,y,cone_origin,x0,y0)

                    
                    isos = vol[:,:,isos_depth].ravel()
                    cost = vol[:,:,cost_depth].ravel()
                    
                    
                    isos_phase = np.angle(isos)
                    cost_phase = np.angle(cost)
                    isos_amp = np.abs(isos)
                    cost_amp = np.abs(cost)
                    isos_real = np.real(isos)
                    cost_real = np.real(cost)
                    isos_imag = np.imag(isos)
                    cost_imag = np.imag(cost)

                    
                    amp = (isos_amp+cost_amp)/2.0

                    valid = np.where(amp>amp.mean()-0*amp.std())[0]
                    d_phase = (cost_phase-isos_phase)
                    d_phase = utils.unwrap(d_phase)
                    d_phase_median = np.median(d_phase[valid])

                    
                    if stimulus:
                        stim_d_phase.append(d_phase_median)
                    else:
                        no_stim_d_phase.append(d_phase_median)

                    
                    #print x0,y0,x,y,fk,ik,stimulus,d_phase_median
            pv = np.nanvar(stim_d_phase)/np.nanvar(no_stim_d_phase)
            print x0,y0,np.var(no_stim_d_phase),np.nanvar(stim_d_phase),pv
            if do_plot and pv>0.01:
                plt.figure(1)
                lu = min(1.0,-np.log(np.var(stim_d_phase)))
                colorVal = scalarMap.to_rgba(lu)
                plt.plot(x0*3,y0*3,marker='o',color=colorVal,alpha=0.5,markersize=10)
                
        if do_plot:
            outfn = '%s_marked_average.png'%self.tag
            plt.figure(1)
            plt.savefig(outfn)
            plt.show()
            
    def crop_average_to_reference(self):
        av = self.h5['/average_image/ISOS_COST'][:,:]
        ref = self.h5['/reference_frame'][:,:]
        x1 = int(self.h5['/reference_coordinates/x1'].value*3)
        x2 = int(self.h5['/reference_coordinates/x2'].value*3)
        y1 = int(self.h5['/reference_coordinates/y1'].value*3)
        y2 = int(self.h5['/reference_coordinates/y2'].value*3)
        av = av[2:,x1-90:x1+510]
        return av


    def make_big_sheet(self,phase=False,sort_by_nvol=True):
        volume_dictionary = self.get_volume_dictionary()
        
        n_volumes = self.get_n_volumes()
        slow_max,fast_max,depth_max = self.get_cone_volume_size()
        n_cones = self.get_n_cones()
        border = 1

        sy = (border+depth_max)*n_cones
        sx = (border+fast_max)*n_volumes
        dpi = 100.0

        fx = sx/dpi
        fy = sy/dpi
            
        cone_catalog = self.h5['cone_catalog']
        cone_keys = cone_catalog.keys()

        clim_dict = {}
        nvol_list = []
        
        for cone_index,ck in enumerate(cone_keys):
            print cone_index
            if ck in ['slow_max','fast_max','depth_max','n_cones']:
                continue
            frame_keys = cone_catalog['%s'%ck].keys()

            cmax = -np.inf
            cmin = np.inf
            n_valid_vols = 0
            
            for fk in frame_keys:
                index_keys = cone_catalog['%s/%s'%(ck,fk)].keys()
                for ik in index_keys:
                    vol = cone_catalog['%s/%s/%s/cone_volume'%(ck,fk,ik)][:,:,:]
                    if np.max(vol.shape)>20:
                        continue
                    if not np.prod(vol.shape):
                        continue
                    if phase:
                        vol = np.angle(vol)
                    else:
                        vol = np.abs(vol)

                    vol = utils.smooth_cone_volume(vol)
                    n_valid_vols = n_valid_vols+1
                    proj = vol.mean(axis=0).T
                    if proj.max()>cmax:
                        cmax = proj.max()
                    if proj.min()<cmin:
                        cmin = proj.min()

            clim_dict[ck] = (cmin,cmax)
            nvol_list.append(n_valid_vols)

        if sort_by_nvol:
            new_cone_keys = []
            order = np.argsort(nvol_list)[::-1]
            for idx in order:
                new_cone_keys.append(cone_keys[idx])
            cone_keys = new_cone_keys
            
            
        sheet = np.zeros((sy,sx),dtype=np.uint8)
        for cone_index,ck in enumerate(cone_keys):
            print cone_index
            if ck in ['slow_max','fast_max','depth_max','n_cones']:
                continue
            frame_keys = cone_catalog['%s'%ck].keys()
            for fk in frame_keys:
                index_keys = cone_catalog['%s/%s'%(ck,fk)].keys()
                for ik in index_keys:
                    vol = cone_catalog['%s/%s/%s/cone_volume'%(ck,fk,ik)][:,:,:]
                    if np.max(vol.shape)>20:
                        continue
                    if not np.prod(vol.shape):
                        continue
                    if phase:
                        vol = np.angle(vol)
                    else:
                        vol = np.abs(vol)
                    vol = utils.smooth_cone_volume(vol)
                    proj = vol.mean(axis=0).T
                    proj = utils.bmpscale(proj,clim_dict[ck])
                    py,px = proj.shape
                    volume_index = volume_dictionary[(fk,ik)]
                    x1 = volume_index*(fast_max+border)
                    x2 = x1+px
                    y1 = cone_index*(depth_max+border)
                    y2 = y1+py
                    sheet[y1:y2,x1:x2] = proj

        plt.figure(figsize=(fx,fy))
        plt.axes([0,0,1,1])
        plt.imshow(sheet,cmap='gray',interpolation='none')
        plt.xticks([])
        plt.yticks([])
        if phase:
            outfn = '%s_cone_catalog_sheet_phase.png'%self.tag
        else:
            outfn = '%s_cone_catalog_sheet_amplitude.png'%self.tag
            
        plt.savefig(outfn,dpi=dpi)
        plt.show()

                            
                
    def get_subvol(self,vol,x,y,z,xrad=0,yrad=0,zrad=0):
        # return a subvol, trimming as necessary
        sx,sy,sz = vol.shape
        x1 = x-xrad
        x2 = x+xrad+1
        y1 = y-yrad
        y2 = y+yrad+1
        z1 = z-zrad
        z2 = z+zrad+1
        x1 = max(0,x1)
        x2 = min(sx,x2)
        y1 = max(0,y1)
        y2 = min(sy,y2)
        z1 = max(0,z1)
        z2 = min(sz,z2)
        return vol[int(x1):int(x2),int(y1):int(y2),int(z1):int(z2)]

    def find_corresponding_data(self,points):
        ref = self.h5['/reference_frame'][:,:]
        #xclicks,yclicks,junk = collector([ref])
        #print xclicks
        #print yclicks

        fkeys = self.h5['/frames'].keys()
        for fk in fkeys:
            ikeys = self.h5['/frames'][fk].keys()
            for ik in ikeys:
                vol = dataset_h5
                print fk,ik
            sys.exit()

        

        
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
        self.h5.put('average_volume/%s_xy_rect'%tag,rect)
        x1 = rect[0]
        x2 = rect[1]
        y1 = rect[2]
        y2 = rect[3]
        
        yproj = np.log(av.mean(axis=0)+1000.0)
        xc,yc,_=collector([yproj])
        zlim = np.round(yc).astype(int)
        self.h5.put('average_volume/%s_zlims'%tag,zlim)
        z1 = zlim[0]
        z2 = zlim[1]
            
        av = av[y1:y2,z1:z2,x1:x2]
        sy,sz,sx = av.shape
        
        bscan = av[sy//2,:,:]
        profile = np.mean(bscan,axis=1)
        lprofile = np.log(profile+1000.0)
        z = np.arange(bscan.shape[0])
        self.h5.put('average_volume/%s_profile'%tag,profile)
        for label in ['ELM','ISOS','ISOS_distal','COST_proximal','COST','RPE','OPL','CH']:
            #xc,zc,_=collector([bscan],titles=['Click upper and lower extents of %s.'%label])
            zc,xc,_=collector([(z,lprofile)],titles=['Click left and right extents of %s.'%label])
            if len(zc)<2:
                sys.exit('Must select two points.')
            zc = np.round(np.array(zc))
            zc = np.array([zc.min(),zc.max()]).astype(int)
            self.h5.put('average_volume/%s_labels/%s'%(tag,label),zc)

    def average_is_cropped_and_labeled(self,tag):
        labeled = True
        try:
            x = self.h5['/average_volume/%s_labels'%tag]
        except Exception as e:
            labeled = False
        return labeled
            
    def get_average_volume_labels(self,tag):
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


    def correct_reference_b(self,do_plot=False,goodness_threshold=0,medfilt_kernel=0,ignore_y=False):

        files = self.h5['frames'].keys()
        frame_keys = []
        for fn in files:
                frames = self.h5['frames/%s'%fn].keys()
                for f in frames:
                    frame_keys.append('/frames/%s/%s'%(fn,f))

        all_yshifts = []
        all_xshifts = []
        all_goodnesses = []
        for fk in frame_keys:
            ys = self.h5[fk]['y_shifts'][:]
            all_yshifts.append(ys)
            xs = self.h5[fk]['x_shifts'][:]
            all_xshifts.append(xs)
            g = self.h5[fk]['goodnesses'][:]
            all_goodnesses.append(g)

        all_yshifts = np.array(all_yshifts)
        all_xshifts = np.array(all_xshifts)
        all_goodnesses = np.array(all_goodnesses)
        # now all_yshifts and all_xshifts are 2D arrays;
        # the first dimension (vertical) corresponds to
        # frame index, and the second dimension (horizontal)
        # corresponds to row index within the frame
        
        if do_plot:
            plt.figure()
            plt.subplot(1,3,1)
            plt.imshow(all_xshifts,interpolation='none',aspect='normal')
            plt.subplot(1,3,2)
            plt.imshow(all_yshifts,interpolation='none',aspect='normal')
            plt.subplot(1,3,3)
            plt.imshow(all_goodnesses,interpolation='none',aspect='normal')

        # diff these:
        # we do this to effectively remove outliers,
        # since the line-to-line differences in x
        # and y can be anticipated while the absolute
        # positions cannot be (as well)
        
        d_yshifts_unfixed = np.diff(all_yshifts,axis=1)
        d_xshifts_unfixed = np.diff(all_xshifts,axis=1)

        # make a goodness mask:
        if goodness_threshold:
            dg = np.min(np.array([all_goodnesses[:,1:],all_goodnesses[:,:-1]]),axis=0)
            goodness_mask = np.ones(dg.shape)*np.nan
            goodness_mask[np.where(dg>goodness_threshold)] = 1.0

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
        # thresholds of .5,-.5 are somewhat artibrary, between
        # thresholds established by maximum drift and saccade speeds
        # Using max values for the latter from Martinez-Conde, 2004
        # (0.5 deg/s and 97 deg/s, resp.), we get limits of
        # 0.01 um and 1.8 um, resp. Practically, any threshold in
        # (0.33,5.0] works, provided we do not exclude
        # the smallest possible shift (0.33 in the case of
        # oversampling by 3)
        def outlier_to_nan(arr,ulim=.34,llim=-.34):
            invalid = np.where(np.logical_or(arr>ulim,arr<llim))
            arr[invalid] = np.nan
            return arr

        d_yshifts = outlier_to_nan(d_yshifts_unfixed)
        d_xshifts = outlier_to_nan(d_xshifts_unfixed)

        if goodness_threshold:
            d_yshifts = d_yshifts*goodness_mask
            d_xshifts = d_xshifts*goodness_mask
        
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
        # original, working version:
        #d_y = np.array([first_y]+list(d_y))
        #d_x = np.array([first_x]+list(d_x))
        # modified version, 2017.06.30
        d_y = np.array([d_y[0]]+list(d_y))
        d_x = np.array([d_x[0]]+list(d_x))


        #if medfilt_kernel:
        #    d_y = medfilt(d_y,medfilt_kernel)
        #    d_x = medfilt(d_x,medfilt_kernel)
        
        if do_plot:
            plt.figure()
            plt.subplot(1,2,1)
            plt.plot(d_x)
            plt.subplot(1,2,2)
            plt.plot(d_y)
            
        vc = np.cumsum(d_y)
        hc = np.cumsum(d_x)
        #if medfilt_kernel:
        #    hc = medfilt(hc,medfilt_kernel)
        #    vc = medfilt(vc,medfilt_kernel)

        if ignore_y:
            vc = vc*0.0

        if do_plot:
            plt.figure()
            plt.subplot(1,3,1)
            plt.plot(hc)
            plt.subplot(1,3,2)
            plt.plot(vc)
            plt.subplot(1,3,3)
            plt.plot(hc,vc)
            

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
                

        canvas_width = int(canvas_width*oversample_factor)
        canvas_height = int((canvas_height+1)*oversample_factor)

        rmean = np.mean(self.reference)

        embedded_reference = np.ones((canvas_height,canvas_width))*rmean
        ref_x1 = int(ref_x1*oversample_factor)
        ref_x2 = int(ref_x2*oversample_factor)
        ref_y1 = int(ref_y1*oversample_factor)
        ref_y2 = int(ref_y2*oversample_factor)
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

        
        canvas_width = int(xmax-xmin+n_fast)
        canvas_height = int(ymax-ymin+n_slow)
        canvas_depth = int(zmax-zmin+max_depth+1)
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

        errorcount = 0
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
                    shifted_bscan[int(z1):int(z2),:] = bscan
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

                try:
                    sum_image[y1:y2,:,x1:x2] = sum_image[y1:y2,:,x1:x2] + block
                    counter_image[y1:y2,:,x1:x2] = counter_image[y1:y2,:,x1:x2] + 1.0
                except Exception as e:
                    errorcount = errorcount + 1
                    print e
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
        print 'error count %d'%errorcount
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
        
            
    def render_stack(self,layer_names):
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

        
        canvas_width = int(xmax-xmin+n_fast)
        canvas_height = int(ymax-ymin+n_slow)
        canvas_depth = int(zmax-zmin+max_depth+1)
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

        errorcount = 0
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
                    shifted_bscan[int(z1):int(z2),:] = bscan
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

                try:
                    sum_image[y1:y2,:,x1:x2] = sum_image[y1:y2,:,x1:x2] + block
                    counter_image[y1:y2,:,x1:x2] = counter_image[y1:y2,:,x1:x2] + 1.0
                except Exception as e:
                    errorcount = errorcount + 1
                    print e
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
        print 'error count %d'%errorcount
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

            
    def filter_registration(self,xshifts,yshifts,goodnesses,xmax=25,ymax=25,medfilt_region=49,do_plot=False):
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

    
