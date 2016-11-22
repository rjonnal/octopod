from octopod import *
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)

fn = '/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.09.20_NFL/15_33_30-JVM_RE_8N_checker.hdf5'

h5 = H5(fn)

ad = ActivityDetector(h5)

ad.project_series()
h5.catalog()

h5.close()
