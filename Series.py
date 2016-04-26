import numpy as np
from matplotlib import pyplot as plt
import sys,os
from octopod import H5
import glob

class Series:

    def __init__(self,reference_h5,vidx=0):

        self.reference_filename = reference_h5
        self.reference_h5 = H5(reference_h5)
        self.reference_vidx = vidx

        self.filenames = [self.reference_filename]
        self.reference_labels = self.reference_h5.get('model')['labels'].keys()


        self.reference_layers = {}

        print self.reference_h5.h5.keys()

        for key in self.reference_labels:
            self.reference_layers[key] = self.reference_h5.get('model/volume_labels/%s'%key)[vidx,:,:]
        
        for key in self.reference_labels:
            plt.figure()
            im = self.reference_layers[key]
            print im.shape
            plt.imshow(im)
            plt.title(key)
        plt.show()

        self.reference_h5.close()


    #def add(self,filename,registration_label):
        
        

    
        



if __name__=='__main__':

    d = '/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.04.12_2/*.hdf5'
    flist = glob.glob(d)
    for f in flist:
        h5 = H5(f)
        print h5.keys()
    
    #s = Series('./oct_test_volume/oct_test_volume_2T.hdf5')
    
