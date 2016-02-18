import sys,os
import h5py
import logging
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from utils import translation,autotrim_bscan,find_peaks
import logging
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
            self.h5 = h5py.File(h5)
            self.logger.info('Opening file %s.'%self.h5)
        else:
            self.h5 = h5
        try:
            self.profile = self.h5['model']['profile'][:]
        except Exception as e:
            self.logger.info('Model does not exist in h5 file.')
            self.profile = self.make_model(debug=debug)

    def blur(self,bscan,kernel_width=5):
        return sp.signal.convolve2d(bscan,np.ones((1,kernel_width)),mode='same')/float(kernel_width)
        
    def make_model(self,vidx=0,debug=False):
        self.logger.info('Making model...')
        avol = np.abs(self.h5['processed_data'][vidx,:,:,:])
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

        model_profile = np.mean(template/counter,axis=1)
        
        self.h5.require_group('model')

        self.write_profile(model_profile)
        
        return model_profile

    def write_profile(self,model_profile):
        try:
            del self.h5['/model/profile']
        except Exception as e:
            pass

        self.h5['model'].create_dataset('profile',data=model_profile)

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
        plt.cla()
        plt.plot(self.profile)
        plt.title('post cropping')
        plt.show()
        
    
    def label(self,smoothing=5):


        # x = np.zeros(10)
        # x[3] = 1.0
        # x[4] = 0.5
        # x[7] = 2.0
        
        # print find_peaks(x,gradient_threshold=1.00001)
        # sys.exit()
        
        if smoothing>1:
            working_profile = sp.signal.convolve(self.profile,np.ones((smoothing)),mode='same')/float(smoothing)
        else:
            working_profile = self.profile
        
        # find peaks and troughs:
        gthresh = 5.0/smoothing

        
        peaks = list(find_peaks(working_profile,gradient_threshold=gthresh))+list(find_peaks(-working_profile,gradient_threshold=gthresh))
        peaks = sorted(peaks)
        
        label_dict = {}

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

        try:
            del self.h5['/model/labels']
        except Exception as e:
            pass

        self.h5['model'].create_group('labels')
        for key in label_dict.keys():
            self.h5['model/labels'].create_dataset(key,data=label_dict[key])

        
