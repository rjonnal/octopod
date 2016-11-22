from octopod import *
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)
import sys,os
import shutil
from glob import glob
import h5py

fm = FileManager()
files = fm.h5_files['2g_aooct']
#files = [f for f in files if f.find('2016.04.12')>-1]

for fn in files:
    m = Model(fn)
    #m.cleanup_modeldb()
    #m.click_label()
    m.write_axial_alignment()