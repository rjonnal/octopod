from octopod import *
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)
import sys,os
import shutil
from glob import glob
import h5py

files = glob('/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.01.29/*.hdf5') + glob('/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.02.02/*.hdf5')

for fn in files:
    h5 = h5py.File(fn)
    p = OCTProcessor(h5)
    p.cleanup()
    p.run()
    h5.close()
    


