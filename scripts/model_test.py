from octopod import Model
import h5py
from matplotlib import pyplot as plt

import glob

#fn = '/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.02.02/3T_vid_1.hdf5'
#fn = '/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.02.02/2T_1.hdf5'

files = glob.glob('/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.01.29/*.hdf5') + glob.glob('/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.02.02/*.hdf5')


for idx,fn in enumerate(files):
    print
    print
    print 'File %d of %d'%(idx+1,len(files))
    model = Model(fn)
    plt.figure()
    plt.plot(model.profile)
    plt.pause(.1)
    
plt.show()
