from octopod import *
from matplotlib import pyplot as plt
import sys,os
import logging
logging.basicConfig(level=logging.INFO)
from octopod.utils import translation
import scipy as sp

fm = FileManager()
files = fm.get('hroct',['Carmen'])

def profile_druse(fn,druse_id):
    dataset = Dataset(fn)
    h5 = dataset.get_h5_handle()
    keys = h5['drusen'][druse_id].keys()
    coordinates = h5['drusen'][druse_id]['coordinates'][:]
    xs = coordinates[:,0]
    ys = coordinates[:,1]
    xmin = int(round(np.min(xs)))
    xmax = int(round(np.max(xs)))
    ymin = int(round(np.min(ys)-150))
    ymax = int(round(np.max(ys)+10))
    reference_index = h5['drusen'][druse_id]['reference_index']
    bscan = np.abs(h5['processed_data'][0,reference_index,:,:])[0,:,:]
    reference = np.abs(h5['processed_data'][0,reference_index,ymin:ymax,xmin:xmax][0,:,:])

#    plt.figure()
#    plt.imshow(bscan)
#    plt.title(os.path.split(fn)[1].replace('.unp',''))
#    plt.autoscale(False)
#    plt.plot([xmin,xmax,xmax,xmin,xmin],[ymax,ymax,ymin,ymin,ymax],'y-')
    
#    plt.figure()
#    plt.imshow(reference,interpolation='none')
#    plt.title(os.path.split(fn)[1].replace('.unp',''))
    sy,sx = reference.shape
    plt.figure()
    smoothed = sp.signal.convolve2d(reference,np.ones((5,1)),'same')
    plt.imshow(smoothed,interpolation='none')
    plt.title(os.path.split(fn)[1].replace('.unp',''))

carmen_ids = ['1_deg_nasal','0p5_deg_temporal','2_deg_temporal','5_deg_temporal']
for f in files:
    profile_druse(f,'2_deg_temporal')

plt.show()


    
#move_coordinates_in_h5('Carmen')
