from octopod import *
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)
import sys,os
import shutil
from glob import glob
import h5py

files = glob('D:/Data/2016.01.29/*.hdf5')+glob('D:/Data/2016.02.02/*.hdf5')+glob('D:/Data/2016.02.24/*.hdf5')

for fn in files:
    source_fn = fn
    path,infn = os.path.split(fn)
    path,datedir = os.path.split(path)
    outroot = 'D:/Data_Share/Dropbox (AO-OCT Review)/Share/2g_aooct_data/Data'
    outpath = os.path.join(outroot,datedir)
    dest_fn = os.path.join(outpath,infn)
    print source_fn,'->',dest_fn
    shutil.move(source_fn,dest_fn)
    
