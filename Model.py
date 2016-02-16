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


    def blur(self,bscan):
        #bprof = np.mean(bscan,axis=1)
        #sy,sx = bscan.shape
        #out = np.tile(bprof,(sx,1)).T
        return sp.signal.convolve2d(bscan,np.ones((1,5)),mode='same')
        
        return out
        
    def make_model(self,vidx=0,debug=False):
        self.logger.info('Making model...')
        avol = np.abs(self.h5['processed_data'][vidx,:,:,:])
        nSlow,nFast,nDepth = avol.shape
        tx_vec = [0.0]
        ty_vec = [0.0]
        g_vec = [1.0]

        if True:
            # make a test volume with very clean translations.
            avol = np.zeros(avol.shape)
            zpos = 100
            rad = 3
            for iSlow in range(nSlow):
                avol[iSlow,zpos-rad:zpos+rad,:] = avol[iSlow,zpos-rad:zpos+rad,:] + 100.0
                zpos = np.abs(np.round(zpos + np.random.randn()))

        template = None

        # compute an intensity threshold for the first template image,
        # to avoid using 

        thresh = np.mean(avol) * .75
        
        for iSlow in range(1,nSlow):
            if debug:
                self.logger.debug('%d of %d'%(iSlow+1,nSlow))

            last = self.blur(avol[iSlow-1,:,:])

            bright_enough = np.mean(last)>=thresh
            
            if template is None:
                if bright_enough:
                    template = last
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
            
            tx,ty,g = translation(target,template)
            print tx,ty,g
            # ty is the amount to shift target to align with template
            # let's say template's profile is [0 1 0 0] and target's is
            # [0 0 1 0]. ty = -1 in this case.

        print tx_vec
        print ty_vec


