import sys,os
import logging
import numpy as np
import scipy as sp
from scipy.ndimage.filters import generic_filter
from scipy.interpolate import bisplrep,bisplev
from matplotlib import pyplot as plt
from scipy.ndimage.morphology import grey_opening
from scipy.ndimage.filters import median_filter
from utils import translation,translation1,autotrim_bscan,find_peaks,shear,Clock,lateral_smooth_3d,polyfit2d,polyval2d
from octopod.Misc import H5
import octopod_config as ocfg
import logging
logging.basicConfig(level=logging.DEBUG)

class Flipper:

    def __init__(self,h5,debug=False):
        """Initialize Flipper. May pass an h5py.File object or .hdf5 filename.
        """
        self.logger = logging.getLogger(__name__)
        if type(h5)==str:
            self.h5 = H5(h5)
            self.logger.info('Opening file %s.'%self.h5)
        else:
            self.h5 = h5

    def flip(self):
        self.logger.info('crop: Starting')
        vols = self.h5.get('/processed_data')[:]

        nv,ns,nd,nf = vols.shape
        subvol = vols[0,:10,:,:]
        test = np.mean(np.abs(subvol),axis=0)
        plt.figure()
        plt.subplot(1,3,1)
        plt.imshow(test)
        plt.subplot(1,3,2)
        plt.imshow(np.log(test))
        plt.subplot(1,3,3)
        plt.plot(np.mean(test,axis=1))
        plt.show()

        flip = raw_input('Flip this set? (y or yes to flip) ')
        flip_me = flip.lower() in ['y','yes']
        if flip_me:
            newvols = np.zeros((nv,ns,nd,nf),dtype=np.complex)
            for v in range(nv):
                print v
                for s in range(ns):
                    old = vols[v,s,:,:]
                    old = np.flipud(old)
                    newvols[v,s,:,:] = old
                    
            subvol = newvols[0,:10,:,:]
            test = np.mean(np.abs(subvol),axis=0)
            
        plt.figure()
        plt.subplot(1,3,1)
        plt.title('post flip; flipped: %s'%flip_me)
        plt.imshow(test)
        plt.subplot(1,3,2)
        plt.imshow(np.log(test))
        plt.title('Ctrl-C now to abort write!')
        plt.subplot(1,3,3)
        plt.plot(np.mean(test,axis=1))
        plt.show()

        if flip_me:
            del self.h5.h5['/processed_data']
            self.h5.put('/processed_data',newvols)
        
