import numpy as np
from octopod import Series, utils
from matplotlib import pyplot as plt
import os,sys


wdir = '/home/rjonnal/Share/2g_aooct_data/Data/2017.08.11/'
fn = os.path.join(wdir,'reg_14_27_33-2.0T_0.0S_no_stimulus_1_000.hdf5')

s = Series(fn)

#print s.h5.catalog(1)
ref = s.h5['reference_frame'][:,:]

cx,cy = utils.find_cones(ref,5,do_plot=False)
cx = np.array(cx)
cy = np.array(cy)
cx = cx - 1
cy = cy - 1

points = zip(cx,cy)

#points = zip([103.88548404780173, 113.65315306336024, 108.07162791161252, 99.30065981600896, 144.15220121398167],[52.054046107384579, 67.801920642672769, 111.05874056871754, 141.75712890332997, 155.91028196669026])

s.find_corresponding_images(points)



