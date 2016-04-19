import glob
import h5py
import numpy as np
from matplotlib import pyplot as plt
import sys,os
from scipy.ndimage.morphology import grey_opening
import scipy as sp


flist = glob.glob('/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.04.12_2/*.hdf5')


for f in flist:
    h5 = h5py.File(f)
    offset_matrix = h5['model/z_offsets'][:]
    goodness_matrix = h5['model/z_offset_goodness'][:]

    om = offset_matrix[0,:,:]
    oms = grey_opening(om,(1,15))


    mode = sp.stats.mode(oms.ravel())[0][0]
    mask = np.zeros(oms.shape)
    lower_threshold = np.mean(oms)-2.0*np.std(oms)
    upper_threshold = np.mean(oms)+2.0*np.std(oms)
    cond = np.logical_and(foffset_matrix>lower_threshold,foffset_matrix<upper_threshold)
    mask[np.where(cond)] = 1


    
    plt.figure(figsize=(18,6))
    clim = np.min(om),np.max(om)
    plt.subplot(131)
    plt.imshow(om,interpolation='none',clim=clim,cmap='gray')
    plt.colorbar()
    plt.subplot(132)
    plt.imshow(oms,interpolation='none',clim=clim,cmap='gray')
    plt.colorbar()
    plt.subplot(133)
    plt.imshow(om-oms,interpolation='none')
    plt.colorbar()


    
plt.show()
    
