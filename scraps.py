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
