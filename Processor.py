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
            

    

