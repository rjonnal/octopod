import numpy as np
from matplotlib import pyplot as plt
import sys,os
from octopod import H5
import glob

class Series:

    def __init__(self,reference_h5_fn,vidx=None,layer_names=['ISOS','COST']):

        self.reference_filename = reference_h5_fn
        self.reference_h5 = H5(reference_h5_fn)
        self.layer_names = layer_names
        
        self.n_vol = self.reference_h5['config']['n_vol'].value
        self.n_fast = self.reference_h5['config']['n_fast'].value
        self.n_slow = self.reference_h5['config']['n_slow'].value

        if vidx is None:
            # if the caller doesn't know which volume to use,
            # show all of them and then quit
            stack = np.zeros((self.n_vol,self.n_slow,self.n_fast))
            for vidx in range(self.n_vol):
                for layer_name in layer_names:
                    stack[vidx,:,:] = stack[vidx,:,:]+self.reference_h5['projections'][layer_name][vidx,:,:]/float(len(layer_names))
            clims = np.percentile(stack[:,50:-50,50:-50],(5,99))
            for vidx in range(self.n_vol):
                plt.figure()
                test = stack[vidx,:,:]
                plt.imshow(test,cmap='gray',interpolation='none',aspect='auto',clim=clims)
                plt.colorbar()
                plt.title('volume %d'%vidx)
            plt.show()
            sys.exit()
                        
        self.reference_vidx = vidx
        
        stack = np.zeros((len(layer_names),self.n_slow,self.n_fast))
        for idx,layer_name in enumerate(self.layer_names):
            stack[idx,:,:] = self.reference_h5['projections'][layer_name][self.reference_vidx,:,:]

        self.reference = np.mean(stack,axis=0)
        plt.imshow(self.reference)
        plt.show()
        self.reference_h5.close()


    #def add(self,filename,registration_label):
        
        

    
        



if __name__=='__main__':

    fn = '/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.08.16/14_44_47-4T_1000.hdf5' # volume 8
    #fn = '/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.08.16/14_44_08-4T_500.hdf5' # volume
    
    s = Series(fn,vidx=8,layer_names=['ISOS'])
    
    #s = Series('./oct_test_volume/oct_test_volume_2T.hdf5')
    
