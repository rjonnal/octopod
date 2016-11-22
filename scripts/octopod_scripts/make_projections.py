from octopod import *
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)
import sys,os
import shutil
from glob import glob
import h5py
import datetime


def now():
    return str(datetime.datetime.now()).replace(' ','_').replace('-','').replace(':','')[:15]
    
logger = logging.getLogger(__name__)
fh = logging.FileHandler('./make_projections_%s.log'%(now()))
fh.setLevel(logging.INFO)
#ch = logging.StreamHandler()
#ch.setLevel(logging.INFO)

#logging.getLogger('').addHandler(ch)
logging.getLogger('').addHandler(fh)

#directories_to_process = ['D:/Data/2015.10.29','D:/Data/2016.01.29','D:/Data/2016.02.02','D:/Data/2015.02.24','D:/Data/2016.03.15-CD','D:/Data/2016.03.15-CW']
#directories_to_process = glob('D:/Data/2016.04.12*')
directories_to_process = glob('D:/Data/*')
#directories_to_process = ['D:/Data/2015.02.24','D:/Data/2016.03.15-CD']+glob('D:/Data/2016.04.12*')


filters_to_dropbox = ['2016.04.12_2/IMG_RE_4TR_12', '2016.04.12_2/IMG_RE_4TR_13', '2016.04.12_2/IMG_RE_4TR_10', '2016.04.12_1/IMG_RE_4TR_1',
                       '2016.04.12_1/IMG_RE_4TR_4', '2016.04.12_1/IMG_RE_4TR_8', '2016.04.12_1/IMG_RE_4TR_10', '2016.04.12_1/IMG_RE_4TR_11',
                       '2016.04.12_1/IMG_RE_4TR_12', '2016.04.12_3']
                       
files = []
for dtp in directories_to_process:
    logger.info('Locating files in directory %s.'%dtp)
    searchstring = os.path.join(dtp,'*.hdf5')
    newfiles = glob(searchstring)
    newfiles = [f.replace('\\','/') for f in newfiles]
    for nf in newfiles:
        logger.info('File to process: %s.'%nf)
    files = files + newfiles


for fidx,hfn in enumerate(files):

    copy_this = False
    for f in filters_to_dropbox:
        if hfn.find(f)>-1:
            copy_this = True
            dest = hfn.replace('D:/Data/','D:/Data_Share/Dropbox/Share/2g_aooct_data/Data/')
    if copy_this:
        destpath,destfilename = os.path.split(dest)
        if not os.path.exists(destpath):
            os.makedirs(destpath)
            print 'making %s'%destpath
    

    logger.info('')
    logger.info('')
    logger.info('Starting work on file %s.'%hfn)
    logger.info('File %d of %d.'%(fidx+1,len(files)))
    h5 = H5(hfn)
    m = Model(h5)
    
    if not h5.has('model/z_offsets') or True:
        m.write_axial_alignment()
        
    m.write_volume_labels()
    
    r = Reporter(h5,hfn.replace('.hdf5','')+'_report')
    try:
        r.projections_report()
    except Exception as e:
        logger.error('Reporter.projections_report EXCEPTION:')
        logger.error('%s'%e)

    if copy_this:
        shutil.copy(hfn,dest)
        logger.info('Copying %s -> %s.'%(hfn,dest))