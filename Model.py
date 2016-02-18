import sys,os
import h5py
import logging
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from utils import translation,autotrim_bscan
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

    def __init__(self,h5):
        """Initialize model. May pass an h5py.File object or .hdf5 filename.
        The file or object must contain a 'processed_data' dataset containing
        at least one volume."""
        if type(h5)==str:
            self.h5 = h5py.File(h5)
        else:
            self.h5 = h5
        self.logger = logging.getLogger(__name__)


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

        plt.figure()
        plt.imshow(template,aspect='auto',interpolation='none')
        plt.colorbar()
        plt.figure()
        plt.imshow(counter,aspect='auto',interpolation='none')
        plt.colorbar()
        plt.show()
