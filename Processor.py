import sys,os
from octopod.DataStore import H5
import logging
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from utils import translation,autotrim_bscan

def process(frame,k_in,k_out,dispersion_coefficients):
    if len(dispersion_coefficients)==2:
        sys.exit('process: did you mean to send 2 dispersion values only?')
    test_frame = frame - np.mean(frame,axis=0)
    test_frame = test_frame.T
    k_interpolator = sp.interpolate.interp1d(k_in,test_frame,axis=0,copy=False)
    test_frame = k_interpolator(k_out)

    dispersion_axis = k_out - np.mean(k_out)
    phase = np.exp(1j*np.polyval(dispersion_coefficients,dispersion_axis))
    test_frame = test_frame * phase[None].T
    test_frame = np.fft.fftshift(np.fft.fft(test_frame,axis=0),axes=0)
    n_depth = test_frame.shape[0]
    test_frame = test_frame[:n_depth/2,:]
    #plt.imshow(np.abs(test_frame))
    #plt.show()
    return test_frame


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
        self.h5.delete(self.post_dataset)


    def run(self):
        # now do something interesting!
        pass

class OCTProcessor(Processor):
    def __init__(self,h5):
        Processor.__init__(self,h5,pre_dataset='/raw_data',post_dataset='/processed_data')
        self.k_in = self.h5.get('k_in')
        self.k_out = self.h5.get('k_out')

    def run(self):
        # if we're re-processing, we should make sure to delete the model and projections
        try:
            self.h5.delete('model')
        except Exception as e:
            self.logger.info('run: Cannot delete model: %s'%e)
        try:
            self.h5.delete('projections')
        except Exception as e:
            self.logger.info('run: Cannot delete projectoins: %s'%e)
            
        n_vol,n_slow,n_fast,n_depth = self.h5.get(self.pre_dataset).shape
        out_block = np.zeros((n_vol,n_slow,n_depth/2,n_fast),dtype=np.complex64)
        c = self.h5.get('dispersion/coefficients')[:]
        for v in range(n_vol):
            for s in range(n_slow):
                print s
                frame = self.h5.get('raw_data')[v][s]
                frame[0,:] = frame[1,:]
                test_frame = process(frame,self.k_in,self.k_out,c)
                out_block[v,s,:,:] = test_frame

        self.h5.put(self.post_dataset,out_block)
        #plt.figure()
        #plt.imshow(np.abs(self.h5.get('processed_data')[0][10]),aspect='auto',interpolation='none')
        #plt.pause(.0001)
        #plt.close()

if __name__=='__main__':

    h5fn = './oct_test_volume/oct_test_volume_2T.hdf5'
    h5 = H5(h5fn)
    op = OCTProcessor(h5)
    op.run()
