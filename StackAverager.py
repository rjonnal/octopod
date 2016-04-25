from octopod import *
import numpy as np
from matplotlib import pyplot as plt
import sys,os,hashlib
import h5py
import logging
logging.basicConfig(level='INFO')

class StackAverager:

    def __init__(self,vol,scale_factor=1.0,cache_dir='./.sa_cache'):
        """Assume the slow dimension is the first dimension."""
        self.vol = vol
        self.id = '%d'%(self.string_to_id(self.vol[0,:,:]))
        self.cache_dir = cache_dir
        if not os.path.exists(cache_dir):
            os.makedirs(cache_dir)
        self.scale_factor = scale_factor
        self.h5 = H5(os.path.join(self.cache_dir,'%s_%0.1fx.hdf5'%(self.id,self.scale_factor)))

    def scale(self,im):
        factor = self.scale_factor
        sy,sx = im.shape
        newshape = int(round(sy*factor)),int(round(sx*factor))
        fim = np.fft.ifft2(np.fft.fftshift(np.fft.fft2(im)),s=newshape)
        aim = np.abs(fim)*factor**2
        return aim

        
    def align_to(self,ref_idx=0,corr_threshold=0.75):
        ref = self.vol[ref_idx,:,:]
        ref = self.scale(ref)
        key = '%03d'%ref_idx

        sy,sx = ref.shape
        try:
            xshifts = self.h5.get('%s/xshifts'%key)[:]
            yshifts = self.h5.get('%s/yshifts'%key)[:]
            corrs = self.h5.get('%s/corrs'%key)[:]
        except Exception as e:
            xshifts = []
            yshifts = []
            corrs = []
            for k in range(self.vol.shape[0]):
                print k
                tar = self.vol[k,:,:]
                tar = self.scale(tar)
                yshift,xshift,corr,nxc = utils.nxcorr2same(tar,ref)
                yshifts.append(yshift)
                xshifts.append(xshift)
                corrs.append(corr)
                
            xshifts = np.array(xshifts).astype(np.int)
            yshifts = np.array(yshifts).astype(np.int)
            corrs = np.array(corrs)

            self.h5.put('%s/xshifts'%key,xshifts)
            self.h5.put('%s/yshifts'%key,yshifts)
            self.h5.put('%s/corrs'%key,corrs)
        

        valid_idx = np.where(corrs>corr_threshold)[0]

        xshifts = xshifts - np.min(xshifts[valid_idx])
        yshifts = yshifts - np.min(yshifts[valid_idx])

        xadd = np.max(xshifts[valid_idx])
        yadd = np.max(yshifts[valid_idx])

        sumimage = np.zeros((ref.shape[0]+yadd,ref.shape[1]+xadd))
        counterimage = np.ones((ref.shape[0]+yadd,ref.shape[1]+xadd))

        for v in valid_idx:
            print v
            xshift = xshifts[v]
            yshift = yshifts[v]
            sumimage[yshift:yshift+sy,xshift:xshift+sx] = sumimage[yshift:yshift+sy,xshift:xshift+sx] + self.vol[v,:,:]
            counterimage[yshift:yshift+sy,xshift:xshift+sx] = counterimage[yshift:yshift+sy,xshift:xshift+sx] + 1

        avimage = sumimage/counterimage
        plt.figure()
        plt.imshow(avimage)
        plt.colorbar()
        plt.figure()
        plt.imshow(sumimage)
        plt.colorbar()
        plt.figure()
        plt.imshow(counterimage)
        plt.colorbar()
        plt.show()
        return
        
        plt.figure()
        plt.plot(corrs)
        
        plt.figure()
        plt.plot(xshifts,'b')
        plt.plot(yshifts,'g')
        plt.show()

    def string_to_id(self,input_string,max_val=2**32):
        md5 = hashlib.md5()
        md5.update(input_string)
        return int(md5.hexdigest(),16)%max_val


    def warp_to_fit(self,ref,tar,strip_width=20):
        sy,sx = ref.shape

        x1s = range(0,sx,strip_width)
        x2s = x1s+strip_width
        x2s[-1] = min(sx,x2s[-1])

        
        
