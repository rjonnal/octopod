import numpy as np
from matplotlib import pyplot as plt
import sys,os
from octopod import H5,utils
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
        self.reference_h5.close()


    def add(self,filename,vidx,slowmin=None,slowmax=None,fastmin=None,fastmax=None):
        if slowmin is None:
            slowmin = 0
        if slowmax is None:
            slowmax = self.n_slow

        if fastmin is None:
            fastmin = 0
        if fastmax is None:
            fastmax = self.n_fast
            
        target_h5 = H5(filename)
        stack = np.zeros((len(self.layer_names),self.n_slow,self.n_fast))
        for idx,layer_name in enumerate(self.layer_names):
            stack[idx,:,:] = target_h5['projections'][layer_name][vidx,:,:]
        target = np.mean(stack,axis=0)
        reference = self.reference
        sy,sx = target.shape
        n = max(sy,sx)
        h = np.hamming(n)
        ham2d = np.sqrt(np.outer(h,h))
        ham2d = ham2d[slowmin:slowmax,fastmin:fastmax]
        y,x,g = utils.strip_register(target[slowmin:slowmax,fastmin:fastmax],reference[slowmin:slowmax,fastmin:fastmax],do_plot=True)
        plt.subplot(2,1,1)
        plt.plot(x)
        plt.plot(y)
        plt.subplot(2,1,2)
        plt.plot(g)
        plt.show()

    
        



if __name__=='__main__':

    fn = '/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.11.21_cones/14_01_03-1T.hdf5' # volume 0
    tfn = '/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.11.21_cones/14_19_25-1T_500_0.hdf5' # volume 0

    
    fn2 = '/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.08.16/14_44_47-4T_1000.hdf5' # volume 8
    #fn = '/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.08.16/14_44_08-4T_500.hdf5' # volume
    
    s = Series(fn,vidx=0,layer_names=['ISOS'])
    s.add(tfn,vidx=5,fastmin=20)
    #s = Series('./oct_test_volume/oct_test_volume_2T.hdf5')
    
