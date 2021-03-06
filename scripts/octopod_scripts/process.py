from octopod import *
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)
import sys,os
import shutil
from glob import glob
import h5py


#directories_to_process = ['D:/Data/2015.10.29','D:/Data/2016.01.29','D:/Data/2016.02.02','D:/Data/2015.02.24','D:/Data/2016.03.15-CD','D:/Data/2016.03.15-CW']
#directories_to_process = glob('D:/Data/2016.04.12*')
#directories_to_process = glob('D:/Data/2016.04.21*')
#directories_to_process = glob('D:/Data/2016.08.12')
#directories_to_process = glob('D:/Data/2016.08.16')
#directories_to_process = glob('D:/Data/2016.09.20_NFL/')
directories_to_process = glob('D:/Data/2016.11.21')

local_path = 'D:/Data/'
cloud_path = 'D:/Data_Share/Dropbox/Share/2g_aooct_data/Data/'

files = []
for dtp in directories_to_process:
    searchstring = os.path.join(dtp,'*.unp')
    print searchstring
    files = files + glob(searchstring)

system_label = '2g_aooct'


# unsupervised steps:
for fn in files:
    continue
    d = Dataset(fn)
    d.initialize(system_label)
    d.optimize_dispersion()
    d.process()

# supervised steps:
for fn in files:
    print fn
    d = Dataset(fn)
    d.crop()
    d.align()
    d.model()

# unsupervised step:
for src in files:
    d = Dataset(src)
    d.label()
    old_path,fn = os.path.split(src)
    new_path = old_path.replace(local_path,cloud_path)
    hsrc = src.replace('.unp','')+'.hdf5'
    junk,hfn = os.path.split(hsrc)
    hdest = os.path.join(new_path,hfn)
    xsrc = src.replace('.unp','')+'.xml'
    junk,xfn = os.path.split(xsrc)
    xdest = os.path.join(new_path,xfn)
    print 'copying %s to %s'%(hsrc,hdest)
    shutil.copyfile(hsrc,hdest)
    print 'copying %s to %s'%(xsrc,xdest)
    shutil.copyfile(xsrc,xdest)
    print
sys.exit()
    
sys.exit()
for fn in files:
    hfn = fn.replace('.unp','.hdf5')
    h5 = H5(hfn)
    r = Reporter(h5,hfn.replace('.hdf5','')+'_report')
    r.processed_report()

sys.exit()    
for fn in files:
    d = Dataset(fn)
    d.label()
    hfn = fn.replace('.unp','.hdf5')
    h5 = H5(hfn)
    r = Reporter(h5,hfn.replace('.hdf5','')+'_report')
    r.processed_report()
    r.dispersion_report()
    r.projections_report(dpi=300)
    