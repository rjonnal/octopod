from octopod import *
import numpy as np
from matplotlib import pyplot as plt
import sys,os,hashlib
import h5py
import logging
logging.basicConfig(level='INFO')

class StackAverager:

    def __init__(self,vol,cache='./.StackAveragerCache.hdf5'):
        """Assume the slow dimension is the first dimension."""

        h5 = h5py.File(cache)
        key = self.string_to_id(vol[0,:,:])

        try:
            self.xshifts = h5['%s/xshifts'][:]
            self.yshifts = h5['%s/yshifts'][:]
            self.corrs = h5['%s/corrs'][:]
        except Exception as e:
            ns,nd,nf = vol.shape
            yshifts = []
            xshifts = []
            corrs = []
            for s in range(1,ns):
                f1 = vol[s-1,:,:]
                f2 = vol[s,:,:]
                yshift,xshift,corr,nxc = utils.nxcorr2same(f1,f2)
                yshifts.append(yshift)
                xshifts.append(xshift)
                corrs.append(corr)

            yshifts = np.array(yshifts).astype(np.int)
            xshifts = np.array(xshifts).astype(np.int)
            corrs = np.array(corrs)
            
            self.yshifts = yshifts
            self.xshifts = xshifts
            self.corrs = corrs
            
            h5.create_dataset('%s/xshifts',data=xshifts)
            h5.create_dataset('%s/yshifts',data=yshifts)
            h5.create_dataset('%s/corrs',data=corrs)
            
        plt.figure()
        plt.plot(corrs)
        
        plt.figure()
        plt.plot(xshifts,'b')
        plt.plot(yshifts,'g')

    def string_to_id(input_string,max_val=2**32):
        md5 = hashlib.md5()
        md5.update(input_string)
        return int(md5.hexdigest(),16)%max_val
