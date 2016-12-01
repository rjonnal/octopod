from octopod import *
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)
import sys,os
import shutil
from glob import glob
import h5py

directories_to_process = glob('/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.11.21_*')

files = []
for dtp in directories_to_process:
    searchstring = os.path.join(dtp,'*.unp')
    files = files + glob(searchstring)

system_label = '2g_aooct'


for idx,f in enumerate(files):
    print idx,f

# unsupervised steps:
for fn in files:
    continue
    d = Dataset(fn)
    d.initialize(system_label)
    #d.optimize_dispersion()
    # hfn = fn.replace('.unp','.hdf5')
    # h5 = H5(hfn)
    # r = Reporter(h5,hfn.replace('.hdf5','')+'_report')
    # r.dispersion_report()
    d.process()
    
# supervised steps:
for fn in files:
    continue
    d = Dataset(fn)
    #d.show()
    #continue
    #d.flip()
    #d.crop()
    #d.align()
    #d.model()

# unsupervised step:
for fn in files:
    d = Dataset(fn)
    d.label()
    hfn = fn.replace('.unp','.hdf5')
    h5 = H5(hfn)
    r = Reporter(h5,hfn.replace('.hdf5','')+'_report')
    r.projections_report(dpi=300)
    
