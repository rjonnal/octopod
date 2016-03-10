from octopod import *
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)
import sys,os
import shutil
from glob import glob
import h5py

files = glob('/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2015.11.17/Cane_Edric/RE/*.unp')

for fn in files:
    if not os.path.exists(fn.replace('.unp','.xml')):
        continue
    ds = Dataset(fn)
    ds.initialize('2g_aooct')
    ds.optimize_dispersion()
    sys.exit()

