import numpy as np
from mayavi import mlab
from octopod import H5,utils
from matplotlib import pyplot as plt
import os,sys


h5 = H5('/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.11.21_cones/14_13_13-1T_500.hdf5')

h5.catalog()

avol = np.abs(h5['flattened_data'][7,:,:,:])
avol = np.transpose(avol,(0,2,1))
mprofile = h5['model']['profile'][:]
avol = avol[10:,25:,:len(mprofile)]

for z in range(avol.shape[2]):
    avol[:,:,z] = utils.gaussian_convolve(avol[:,:,z],0.25)

profile = np.mean(np.mean(avol,axis=1),axis=0)
offset = utils.translation1(profile,mprofile)[0]

costidx = h5['model']['labels']['COST'].value + offset
isosidx = h5['model']['labels']['ISOS'].value + offset

rad = 3

isos = np.mean(avol[:,:,isosidx-rad:isosidx+rad+1],axis=2)
cost = np.mean(avol[:,:,costidx-rad:costidx+rad+1],axis=2)
cones = np.mean(avol[:,:,isosidx-rad:costidx+rad+1],axis=2)

def remove_overlaps(cx,cy,thresh=2):
    N = len(cx)
    cxo,cyo = [],[]
    
    for a in range(N-1):
        safe = True
        for b in range(a+1,N):
            ya = cy[a]
            yb = cy[b]
            xa = cx[a]
            xb = cx[b]
            d = np.sqrt((ya-yb)**2+(xa-xb)**2)
            if d<=thresh:
                safe = False
                break
        if safe:
            cxo.append(cx[a])
            cyo.append(cy[a])

    return np.array(cxo),np.array(cyo)


cones = utils.background_subtract(cones)

# find the cones
nstd = 0.25
left = cones[1:-1,1:-1] - cones[1:-1,:-2]
right = cones[1:-1,2:] - cones[1:-1,1:-1]
top = cones[1:-1,1:-1] - cones[:-2,1:-1]
bottom = cones[2:,1:-1] - cones[1:-1,1:-1]

cones_center = cones[1:-1,1:-1]

htest = np.logical_and(left>0,right<0)
vtest = np.logical_and(top>0,bottom<0)

peak_test = np.logical_and(htest,vtest)
cones_bright_test = cones_center>np.median(cones_center)+nstd*np.std(cones_center)
test = peak_test*cones_bright_test

cy,cx = np.where(test)

cy = cy + 1
cx = cx + 1

cx,cy = remove_overlaps(cx,cy)

newvol = np.zeros(avol.shape)

offsets = [-1,0,1]
fracs = [0.5,1.0,0.5]

for x,y in zip(cx,cy):
    # get the profile:
    prof = avol[y,x,:]
    this_isos_idx = utils.ascend(prof,isosidx+1)
    this_cost_idx = utils.ascend(prof,costidx)

    for o,f in zip(offsets,fracs):
        newvol[y,x,this_isos_idx+o] = f*prof[this_isos_idx]
        newvol[y,x,this_cost_idx+o] = f*prof[this_cost_idx]


    
    # plt.cla()
    # plt.plot(prof)
    # plt.plot(this_isos_idx,prof[this_isos_idx],'ks')
    # plt.plot(this_cost_idx,prof[this_cost_idx],'ks')
    # plt.pause(1)

for z in range(newvol.shape[2]):
    newvol[:,:,z] = utils.gaussian_convolve(newvol[:,:,z],0.5)

#newvol = newvol + avol/2.0
    
nisos = np.mean(newvol[:,:,isosidx-rad:isosidx+rad+1],axis=2)
ncost = np.mean(newvol[:,:,costidx-rad:costidx+rad+1],axis=2)
ncones = np.mean(newvol[:,:,isosidx-rad:costidx+rad+1],axis=2)

    
plt.imshow(ncones,cmap='gray')#,interpolation='none')
plt.autoscale(False)
plt.plot(cx,cy,'b.')
plt.show()

scatter_field = mlab.pipeline.scalar_field(newvol)
vmin,vmax = np.percentile(ncones,[5,99])
mlab.pipeline.volume(scatter_field,vmin=vmin,vmax=vmax)
mlab.pipeline.image_plane_widget(scatter_field,
                                 plane_orientation='z_axes',
                                 slice_index=avol.shape[2]-costidx,
)

# mlab.pipeline.iso_surface(scatter_field, contours=[vmax, ],)
# mlab.pipeline.image_plane_widget(scatter_field,
#                     plane_orientation='z_axes',
#                     slice_index=10)

mlab.show()

sys.exit()



mlab.figure(1, bgcolor=(0, 0, 0), size=(350, 350))
mlab.clf()

# The position of the atoms
atoms_x = np.array([2.9, 2.9, 3.8]) * 40 / 5.5
atoms_y = np.array([3.0, 3.0, 3.0]) * 40 / 5.5
atoms_z = np.array([3.8, 2.9, 2.7]) * 40 / 5.5

O = mlab.points3d(atoms_x[1:-1], atoms_y[1:-1], atoms_z[1:-1],
                  scale_factor=3,
                  resolution=20,
                  color=(1, 0, 0),
                  scale_mode='none')

H1 = mlab.points3d(atoms_x[:1], atoms_y[:1], atoms_z[:1],
                   scale_factor=2,
                   resolution=20,
                   color=(1, 1, 1),
                   scale_mode='none')

H2 = mlab.points3d(atoms_x[-1:], atoms_y[-1:], atoms_z[-1:],
                   scale_factor=2,
                   resolution=20,
                   color=(1, 1, 1),
                   scale_mode='none')

# The bounds between the atoms, we use the scalar information to give
# color
mlab.plot3d(atoms_x, atoms_y, atoms_z, [1, 2, 1],
            tube_radius=0.4, colormap='Reds')

# Display the electron localization function ##################################

# Load the data, we need to remove the first 8 lines and the '\n'
str = ' '.join(file('h2o-elf.cube').readlines()[9:])
data = np.fromstring(str, sep=' ')
data.shape = (40, 40, 40)

source = mlab.pipeline.scalar_field(data)
min = data.min()
max = data.max()
vol = mlab.pipeline.volume(source, vmin=min + 0.65 * (max - min),
                                   vmax=min + 0.9 * (max - min))

mlab.view(132, 54, 45, [21, 20, 21.5])

mlab.show()
