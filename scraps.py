    def render_stack(self,layer_names):
        files = self.hive['frames'].keys()
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
            keys = self.hive['/frames/%s'%filename].keys()
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
                
                goodnesses = self.hive['/frames/%s/%s/goodnesses'%(filename,k)][:]
                xshifts = sign*self.hive['/frames/%s/%s/x_shifts'%(filename,k)][:]
                yshifts = sign*self.hive['/frames/%s/%s/y_shifts'%(filename,k)][:]

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
        
        self.hive.put('/reference_coordinates/x1',ref_x1)
        self.hive.put('/reference_coordinates/x2',ref_x2)
        self.hive.put('/reference_coordinates/y1',ref_y1)
        self.hive.put('/reference_coordinates/y2',ref_y2)
        
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
        self.hive.put('/counter_volume/%s'%label,counter_image)
        self.hive.put('/average_volume/%s'%label,av)
        self.hive.put('/sum_volume/%s'%label,sum_image)
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

            
###################################################################################
            continue

                        try:
                            fit,fit_covar = c.gaussian_mixture_model(True)
                        except Exception as e:
                            print e
                            fit,fit_covar = c.gaussian_mixture_model(True)
                            
                        fit_isos_z,_,fit_cost_z,fit_isos_sigma,_,fit_cost_sigma,fit_isos_amplitude,_,fit_cost_amplitude = fit
                        
                        perr = np.sqrt(np.diag(fit_covar))
                        fit_isos_z_error,_,fit_cost_z_error,fit_isos_sigma_error,_,fit_cost_sigma_error,fit_isos_amplitude_error,_,fit_cost_amplitude_error = perr







def find_corresponding_images_old(self,points,minimum_goodness=10.0,output_radius=2,match_radius=2.0,do_plot=False):
        
        fkeys = self.h5['/frames'].keys()
        for fk in fkeys:
            dataset_fn = os.path.join(self.working_directory,fk)
            dataset_h5 = H5(dataset_fn)
            ikeys = self.h5['/frames'][fk].keys()
            for ik in ikeys:
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

                if do_plot:
                    plt.subplot(1,2,1)
                    plt.imshow(self.reference,cmap='gray',interpolation='none')
                # # quickly project the target volume
                cmed = int(np.median(cost_depths))
                imed = int(np.median(isos_depths))
                # aproj = np.mean(np.mean(np.abs(vol),axis=2),axis=0)
                # plt.plot(aproj)
                # plt.axvline(imed)
                # plt.axvline(cmed)
                # plt.title('%d,%d'%(imed,cmed))
                # plt.show()
                # continue
                proj = np.abs(vol[:,imed-2:cmed+2,:]).mean(axis=1)

                if do_plot:
                    plt.subplot(1,2,2)
                    plt.imshow(proj,cmap='gray',interpolation='none')
                    colors = 'rgbcymkw'

                plt.figure(figsize=(18,9))
                for idx,(ptx,pty) in enumerate(points):
                    if do_plot:
                        color = colors[idx%len(colors)]
                        plt.subplot(1,2,1)
                        plt.plot(ptx,pty,'%so'%color)
                    yerr = np.abs(pty-yramp)
                    match_index = np.argmin(yerr)
                    if do_plot:
                        plt.title(g[match_index])
                    if yerr[match_index]<=match_radius and g[match_index]>minimum_goodness:
                        yout = int(match_index)
                        xout = int(ptx + x[match_index])
                        xout,yout = utils.ascend2d(proj,xout,yout,do_plot=False)
                        
                        cone_volume = self.get_subvol(vol,yout,0,xout,output_radius,np.inf,output_radius)

                        cone_profile = np.mean(np.mean(np.abs(cone_volume),axis=2),axis=0)

                        # let's try to label this cone's peaks
                        shift,corr = utils.nxcorr(model_profile,cone_profile)
                        isos_guess = int(model_isos + shift)
                        cost_guess = int(model_cost + shift)

                        z1 = isos_guess-3
                        z2 = cost_guess+3
                        
                        sy,sz,sx = cone_volume.shape

                        aline_corrs = []
                        dpeaks = []
                        
                        for idx,coney in enumerate(range(sy)):
                            
                            max_displacement = 4
                            cut_profile = np.abs(cone_volume[coney,:,:]).mean(axis=1)

                            # how representitive is it?
                            aline_corr = np.corrcoef(np.vstack((cut_profile[z1:z2+1],cone_profile[z1:z2+1])))[0,1]

                            peaks = utils.find_peaks(cut_profile)
                            heights = cut_profile[peaks]/cut_profile.std()
                            isos_dpeaks = np.abs(peaks-isos_guess)
                            isos_dpeaks[np.where(isos_dpeaks>=max_displacement)] = 2**16
                            isos_scores = heights - isos_dpeaks

                            cost_dpeaks = np.abs(peaks-cost_guess)
                            cost_dpeaks[np.where(cost_dpeaks>=max_displacement)] = 2**16
                            cost_scores = heights - cost_dpeaks

                            isos_peak = peaks[np.argmax(isos_scores)]
                            cost_peak = peaks[np.argmax(cost_scores)]

                            dpeaks.append(cost_peak-isos_peak)
                            aline_corrs.append(aline_corr)
                            
                            if do_plot:
                                plt.subplot(5,3,idx*3+1)
                                plt.imshow(np.mean(np.abs(cone_volume),axis=0),cmap='gray',aspect='normal')
                                plt.ylim((cost_guess+4,isos_guess-4))
                                #plt.ylim((cost_guess-5,isos_guess+5))
                                plt.subplot(5,3,idx*3+2)
                                plt.imshow(np.mean(np.abs(cone_volume),axis=2).T,cmap='gray',aspect='normal')
                                plt.ylim((cost_guess+4,isos_guess-4))
                                #plt.ylim((cost_guess-5,isos_guess+5))

                                plt.subplot(5,3,idx*3+3)
                                plt.plot(cut_profile)
                                plt.plot(cone_profile)
                                plt.axvline(isos_guess,color='g')
                                plt.axvline(cost_guess,color='g')
                                plt.axvline(isos_peak,color='r')
                                plt.axvline(cost_peak,color='r')
                                plt.xlim((isos_guess-4,cost_guess+4))
                                #print cost_peak-isos_peak,aline_corr

                        # computed weighted average of dpeaks:
                        aline_corrs = np.array(aline_corrs)
                        dz = int(np.sum(dpeaks*aline_corrs)/np.sum(aline_corrs))

                        #probe = np.zeros(z2-z1+1)
                        #probe[0] = 1.0
                        #probe[weighted_average] = 1.0

                        def mysum(x,y):
                            return x+y
                            #return (x+y)/2.0
                            #return x+y+np.sqrt(x*y)
                        colors = 'rgbkcym'
                        plt.cla()
                        sheet = []
                        unaligned_sheet = []
                        for idx,coney in enumerate(range(sy)):
                            color = colors[idx%len(colors)]
                            cut_profile = np.abs(cone_volume[coney,:,:]).mean(axis=1)
                            #shift,corr = utils.nxcorr(cut_profile[z1:z2+1],probe)
                            totals = []
                            zrange = np.arange(z1,z1+6)
                            for k in zrange:
                                totals.append(mysum(cut_profile[k],cut_profile[k+dz]))

                            plt.plot(zrange,totals,'%s-'%color)
                            
                            isos = z1+np.argmax(totals)
                            cost = isos+dz
                            cut = cone_volume[coney,isos-3:cost+4,:]
                            sheet.append(cut)

                            unaligned_cut = cone_volume[coney,z1:z2,:]
                            unaligned_sheet.append(unaligned_cut)
                            

                        sheet = np.hstack(sheet)
                        unaligned_sheet = np.hstack(unaligned_sheet)

                        uamplitude = np.abs(unaligned_sheet)
                        
                        amplitude = np.abs(sheet)
                        phase = np.angle(sheet)

                        for row in [4,-4]:
                            phase[row] = unwrap(phase[row])
                        
                        mamp = amplitude.mean(axis=0)
                        valid = np.where(mamp>mamp.mean())[0]
                        
                        isos_phase = phase[3,:]
                        phase = (phase-isos_phase)%(2*np.pi)
                        phase_std = phase[:,valid].std(axis=1)

                        plt.subplot(2,2,1)
                        plt.imshow(amplitude,cmap='gray')
                        plt.colorbar()

                        plt.subplot(2,2,2)
                        plt.imshow(uamplitude,cmap='gray')
                        plt.colorbar()
                        
                        plt.subplot(2,2,3)
                        plt.imshow(phase,cmap='jet')
                        plt.colorbar()
                        
                        plt.subplot(2,2,4)
                        plt.plot(phase_std)
                        plt.show()





                        # scan through the volume and check correlation w/ en face projection to identify
                        # most likely cone layers
                        efp = np.median(np.abs(output_volume[:,isos_guess-2:cost_guess+2,:]),axis=1)
                        corr_vec = []
                        for zpos in range(isos_guess-2,cost_guess+2+1):
                            test = np.abs(output_volume[:,zpos,:]).ravel()
                            corr_vec.append(np.corrcoef(np.vstack((efp.ravel(),test)))[0,1])

                        plt.plot(corr_vec)
                        plt.show()
                        continue
                        
########### Code from Series.py #######################################################################

    
    def render_old(self,layer_names=None,goodness_threshold=0.0,correlation_threshold=0.0,overwrite=False,oversample_factor=3,do_plot=False,frames_to_save=[]):

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

        for filename in files:
            keys = self.h5['/frames/%s'%filename].keys()
            for k in keys:
                test = self.get_image(filename,0,None)
                n_slow,n_fast = test.shape

                goodnesses = self.h5['/frames/%s/%s/goodnesses'%(filename,k)][:]
                xshifts = sign*self.h5['/frames/%s/%s/x_shifts'%(filename,k)][:]
                yshifts = sign*self.h5['/frames/%s/%s/y_shifts'%(filename,k)][:]

                xshifts = np.squeeze(xshifts)
                yshifts = np.squeeze(yshifts)

                xshifts,yshifts,goodnesses = self.filter_registration_2(xshifts,yshifts,goodnesses)

                yshifts = yshifts + np.arange(n_slow)

                assert len(yshifts)==n_slow
                
                valid = np.where(goodnesses>=goodness_threshold)[0]

                print len(valid)
                
                if len(valid):
                    t = np.arange(len(xshifts))
                    xshifts = xshifts[valid]
                    yshifts = yshifts[valid]
                    t = t[valid]

                    #plt.figure()
                    #plt.plot(t,xshifts)
                    #plt.plot(t,yshifts-t)
                    #plt.show()
                
                    newxmin = np.min(xshifts)
                    newxmax = np.max(xshifts)
                    newymin = np.min(yshifts)
                    newymax = np.max(yshifts)

                    xmin = min(xmin,newxmin)
                    xmax = max(xmax,newxmax)
                    ymin = min(ymin,newymin)
                    ymax = max(ymax,newymax)



        print ymin,ymax
        sys.exit()
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

        print embedded_reference.shape,ery1,ery2
        print ref_oversampled.shape
        sys.exit()
        
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
                    current_sum_image[y1:y2,x1:x2] = current_sum_image[y1:y2,x1:x2] + block
                    current_counter_image[y1:y2,x1:x2] = current_counter_image[y1:y2,x1:x2] + 1.0
                    correlation_vector = np.ones(goodnesses.shape)
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
                            x1 = int(np.round(xshifts[v]*oversample_factor-xoffset))
                            y1 = int(v*oversample_factor+np.round(yshifts[v]*oversample_factor-yoffset))
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

                            ref_section = embedded_reference[y1:y2,x1:x2].ravel()
                            corr_valid = np.where(ref_section)[0]
                            a,b = ref_section[corr_valid],block.ravel()[corr_valid]
                            corr = np.corrcoef(a,b)[1,0]

                            correlation_image[y1:y2,x1:x2] = correlation_image[y1:y2,x1:x2] + corr
                            correlation_vector[v] = corr
                            if corr>correlation_threshold:
                                sum_image[y1:y2,x1:x2] = sum_image[y1:y2,x1:x2] + block
                                counter_image[y1:y2,x1:x2] = counter_image[y1:y2,x1:x2] + 1.0
                                current_sum_image[y1:y2,x1:x2] = current_sum_image[y1:y2,x1:x2] + block
                                current_counter_image[y1:y2,x1:x2] = current_counter_image[y1:y2,x1:x2] + 1.0

                self.h5.put('/frames/%s/%s/correlations'%(filename,k),correlation_vector)

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
            
        self.h5.put('%s/correlation_image'%corrstr,correlation_image)
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
        

    def filter_registration_old(self,xshifts,yshifts,goodnesses,xmax=10,ymax=10):

        xvalid = np.abs(xshifts - np.median(xshifts))<=xmax
        yvalid = np.abs(yshifts - np.median(yshifts))<=ymax
        invalid = np.logical_or(1-xvalid,1-yvalid)
        goodnesses[invalid] = -np.inf

        return xshifts,yshifts,goodnesses
        
    def filter_registration_2(self,xshifts,yshifts,goodnesses,xmax=3,ymax=3):

        xmed = medfilt(xshifts,9)
        ymed = medfilt(yshifts,9)
        xerr = np.abs(xshifts-xmed)
        yerr = np.abs(yshifts-ymed)


        plt.subplot(1,2,1)
        plt.plot(xshifts)
        plt.plot(yshifts)
        plt.plot(xmed)
        plt.plot(ymed)
        plt.subplot(1,2,2)
        plt.plot(xerr)
        plt.plot(yerr)
        plt.show()
        
        xvalid = xerr<=xmax
        yvalid = yerr<=ymax

        invalid = np.logical_or(1-xvalid,1-yvalid)
        goodnesses[invalid] = -np.inf

        return xshifts,yshifts,goodnesses
        


########### Code from original StripRegistrar #########################################################

            #nxcval = np.real(np.fft.fftshift(np.fft.ifft2(np.fft.fft2(strip,s=(sy,sx))*np.conj(np.fft.fft2(self.ref)))))
            #sf = np.fft.fftshift(np.fft.fft2(strip,s=(sy,sx)))
            #rf = np.fft.fftshift(np.conj(np.fft.fft2(self.ref)))
            #nxcval = np.fft.ifft2(sf*rf)
            #nxcval1 = np.abs(nxcval)

            # use built-in fft-based convolution, but have to flip ref horizontally and vertically for
            # correlation
            print yshift,xshift,corr
            plt.cla()
            plt.imshow(nxcval,interpolation='none',cmap='gray')
            plt.pause(.1)
            continue
            # #nxcval = utils.background_subtract(nxcval,utils.strel(diameter=25))
            # #nxcval = correlate2d(strip,ref)
            

            # peakVal = np.max(nxcval)
            # peakIdx = np.where(nxcval==peakVal)

            # yoff,xoff = peakIdx[0][0],peakIdx[1][0]


            # print yoff,xoff
            
            # if nxcval.shape[0]%2:
            #     yshift = (nxcval.shape[0]-1)/2.0 - yoff
            # else:
            #     yshift = nxcval.shape[0]/2.0 - yoff

            # yshift = yshift%sy-y1
            

            # if nxcval.shape[1]%2:
            #     xshift = (nxcval.shape[1]-1)/2.0 - xoff
            # else:
            #     xshift = nxcval.shape[1]/2.0 - xoff

            #corrs.append(peakVal/ref.shape[1]/ref.shape[0]/self.strip_width)



            # from when i was trying to insert the reference image in the stack:
        yoffset = -np.min(y1s)
        xoffset = -np.min(x1s)

        if yoffset>=0:
            # all strips are below top of reference
            ry1 = 0
            ry2 = self.rsy
            outsy = max(np.max(y2s)+yoffset,ry2)
        if xoffset>=0:
            # all strips are to the right of reference
            rx1 = 0
            rx2 = self.rsx
            outsx = max(np.max(x2s)+xoffset,rx2)
        if yoffset<0:
            # some strips are above top of reference
            ry1 = -yoffset
            ry2 = ry1 + self.rsy
            outsy = max(np.max(y2s)+yoffset,ry2)
        if xoffset<0:
            # some strips start to the left of reference left edge
            rx1 = -xoffset
            rx2 = rx1 + self.rsx
            outsx = max(np.max(x2s)+xoffset,rx2)

            

            
########### Code from original Model.py in which median filtering was used to smooth the offset matrix
        # first make a 
        




        
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



def crop(self):
        plt.plot(self.profile)
        plt.title('enter cropping coordinates in console')
        plt.pause(.1)

        z1default = np.where(self.profile)[0][0]
        z2default = np.where(self.profile)[0][-1]
        
        z1str = raw_input('Enter starting index of valid region [%d]: '%z1default)
        z2str = raw_input('Enter ending index of valid region: [%d]:'%z2default)
        
        plt.close()
        
        if not len(z1str):
            z1 = z1default
        else:
            z1 = int(z1str)

        if not len(z2str):
            z2 = z2default
        else:
            z2 = int(z2str)
        
        self.profile[:z1] = 0.0
        self.profile[z2+1:] = 0.0
        self.write_profile(self.profile)
        plt.close()
        
    
    def label(self,smoothing=5):
        if smoothing>1:
            working_profile = sp.signal.convolve(self.profile,np.ones((smoothing)),mode='same')/float(smoothing)
        else:
            working_profile = self.profile
        
        # find peaks and troughs:
        gthresh = 1.0/smoothing

        peaks = list(find_peaks(working_profile,gradient_threshold=gthresh))+list(find_peaks(-working_profile,gradient_threshold=gthresh))
        peaks = sorted(peaks)
        
        label_dict = self.get_label_dict()
        
        idx = 0
        z = np.arange(len(working_profile))
        done = False or not len(peaks)
        
        while not done:
            
            peak = peaks[idx]
            
            z1 = max(0,peak-10)
            z2 = min(len(working_profile),peak+10)
            
            plt.subplot(1,2,1)
            plt.cla()
            plt.plot(working_profile)
            plt.plot(peak,working_profile[peak],'k*')
            plt.plot(self.profile)
            plt.plot(peak,self.profile[peak],'go')
            plt.subplot(1,2,2)
            plt.cla()
            plt.plot(z[z1:z2],working_profile[z1:z2])
            plt.plot(peak,working_profile[peak],'k*')
            plt.plot(z[z1:z2],self.profile[z1:z2])
            plt.plot(peak,self.profile[peak],'go')

            
            plt.pause(.1)
            label = raw_input('Enter label for marked peak [q for quit]: ')
            
            if label=='':
                idx = (idx + 1)%len(peaks)
            elif label[0]=='+' or label[0]=='-':
                idx = (idx + int(label))%len(peaks)
            elif label=='q':
                done = True
            else:
                label_dict[label] = peak

        plt.close()

        try:
            del self.h5['/model/labels']
        except Exception as e:
            pass

        self.h5['model'].create_group('labels')
        for key in label_dict.keys():
            self.h5['model/labels'].create_dataset(key,data=label_dict[key])

            prof = np.mean(im2,axis=1)
            newprof = np.zeros(prof.shape)
            offset = (-1)*ty
            print offset
            if offset<0:
                newprof[-offset:] = prof[:offset]
                if False:
                    plt.figure(2)
                    plt.plot(self.prof/count)
                    plt.plot(newprof)
                    plt.pause(.5)
            elif offset>0:
                newprof[:-offset] = prof[offset:]
                if False:
                    plt.figure(2)
                    plt.plot(self.prof/count)
                    plt.plot(newprof)
                    plt.show(.5)
            else:
                newprof = prof
            self.prof = self.prof + newprof
            count = count + 1
            
            if debug:
                plt.figure(1)
                plt.cla()
                plt.plot(self.prof)
                #plt.imshow(im1,aspect='auto',interpolation='none')
                plt.pause(.0001)
