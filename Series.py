import numpy as np
from matplotlib import pyplot as plt
import sys,os
from octopod import H5

class Series:

    def __init__(self,reference_h5):
        if type(reference_h5)==str:
            self.reference_h5 = H5(reference_h5)
        else:
            self.reference_h5 = reference_h5

        self.labels = self.reference_h5.get('model')['labels'].keys()

        print self.labels
        



if __name__=='__main__':

    s = Series('./oct_test_volume/oct_test_volume_2T.hdf5')
