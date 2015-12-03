import sys,os
import h5py
import logging
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt

class Processor:
    """An interface for a data processing class. Classes which
    subclass this should be thought of as processing plugins.
    Some rules:
    1. The Processor's run method reads data from pre_dataset and writes
       data to post_dataset.
    2. The cleanup method deletes the post_dataset, if it
       exists.
    """
    
    def __init__(self,h5,pre_dataset='/raw_data',post_dataset='/test'):
        self.h5 = h5
        self.pre_dataset = pre_dataset
        self.post_dataset = post_dataset
        self.cleanup()

    def cleanup(self):
        try:
            del self.h5[self.post_dataset]
        except Exception as e:
            pass


    def run(self):
        # now do something interesting!
        pass

class ProcessorCopyRaw(Processor):
    """A toy subclass of Processor."""
    def __init__(self,h5):
        Processor.__init__(self,h5,pre_dataset='/raw_data',post_dataset='/copy_of_raw')
        #self.pre_dataset = pre_dataset
        #self.post_dataset = post_dataset
        
    def run(self):
        n_vol,n_slow,n_fast,n_depth = self.h5[self.pre_dataset].shape

        self.h5.create_dataset(self.post_dataset,(n_vol,n_slow,n_fast,n_depth))
        for i_vol in range(n_vol):
            self.h5[self.post_dataset][i_vol] = self.h5[self.pre_dataset][i_vol]
            


class OCTProcessor(Processor):
    def __init__(self,h5):
        Processor.__init__(self,h5,pre_dataset='/raw_data',post_dataset='/processed_data')
        self.k_in = self.h5['k_in']
        self.k_out = self.h5['k_out']

    def run(self):
        n_vol,n_slow,n_fast,n_depth = self.h5[self.pre_dataset].shape
        out_block = np.zeros((n_vol,n_slow,n_depth/2,n_fast),dtype=np.complex64)
        c = self.h5['dispersion/coefficients'][:]
        for v in range(n_vol):
            for s in range(n_slow):
                print s
                frame = self.h5['raw_data'][v][s]
                test_frame = frame - np.mean(frame,axis=0)
                test_frame = test_frame.T
                k_interpolator = sp.interpolate.interp1d(self.k_in,test_frame,axis=0,copy=False)
                test_frame = k_interpolator(self.k_out)

                dispersion_axis = self.k_out - np.mean(self.k_out)
                phase = np.exp(1j*np.polyval(c,dispersion_axis))
                test_frame = test_frame * phase[None].T
                test_frame = np.fft.fftshift(np.fft.fft(test_frame,axis=0),axes=0)
                test_frame = test_frame[:n_depth/2,:]
                out_block[v,s,:,:] = test_frame

        self.h5.create_dataset(self.post_dataset,data=out_block,dtype='c8')
        plt.imshow(np.abs(self.h5['processed_data'][0][10]))
        plt.show()
