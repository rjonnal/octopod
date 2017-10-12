import numpy as np
from octopod import Series, utils
from matplotlib import pyplot as plt
import os,sys


wdir = '/home/rjonnal/Share/2g_aooct_data/Data/2017.08.09/2.0T_0.0S_f0'
hive_name = os.path.join(wdir,'reg_14_27_33-2.0T_0.0S_no_stimulus_1_000')
#fn = os.path.join(wdir,'reg_stim_with_other_2.0T_0.0S_f0_000.hdf5')

s = Series(hive_name)

try:
    fudge
    cone_catalog = s.hive['cone_catalog']
except Exception as e:
    print e
    ref = s.hive['reference_frame']
    print ref
    
    cx,cy = utils.find_cones(ref,5,do_plot=False)
    cx = np.array(cx)
    cy = np.array(cy)
    cx = cx - 1
    cy = cy - 1
    plt.imshow(ref,cmap='gray',interpolation='none')
    plt.autoscale(False)
    plt.plot(cx,cy,'r+')
    plt.show()
    points = zip(cx,cy)
    s.make_cone_catalog(points)

#s.make_big_sheet(phase=False)
#s.analyze_cone_phase()
#s.crop_average_to_reference()
sys.exit()

cone_keys = cone_catalog.keys()
max_length = 0

n_volumes = s.count_volumes()
print n_volumes
sys.exit()

for ck in cone_keys:
    xy = ck.split('_')
    refx = int(xy[0])
    refy = int(xy[1])
    frame_keys = cone_catalog['%s'%ck].keys()
    for fk in frame_keys:
        index_keys = cone_catalog['%s/%s'%(ck,fk)].keys()
        for ik in index_keys:
            dims = cone_catalog['%s/%s/%s/cone_volume'%(ck,fk,ik)].shape
            print ck,fk,ik,dims
            if dims[2]>max_length:
                max_length = dims[2]
                
sys.exit()

#print s.h5.catalog(1)



