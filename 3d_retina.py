import numpy as np
from mayavi import mlab
from octopod import H5,utils
from matplotlib import pyplot as plt
import os,sys
import cones
from cone_density import ConeDensityInterpolator
import octopod_config as ocfg

h5 = H5('/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.11.21_cones/14_13_13-1T_500.hdf5')
vidx = 7

M_PER_DEG = 300e-6

def get_Sz(L1,L2,n=1.38):
    return np.abs(1.0/((2.0*n)/L1 - (2.0*n)/L2))

spectrum_polynomial = ocfg.source_spectra['2g_aooct']
spectrum = np.polyval(spectrum_polynomial,np.arange(736))
Sz = get_Sz(spectrum[0],spectrum[-1])

cdi = ConeDensityInterpolator()

n_vol = h5['config']['n_vol'].value
n_slow = h5['config']['n_slow'].value
n_fast = h5['config']['n_fast'].value
n_depth = h5['config']['n_depth'].value

nt = -h5['eccentricity']['nasal_temporal'].value
si = h5['eccentricity']['superior_inferior'].value

x_scan_mv = h5['config']['x_scan_mv'].value
y_scan_mv = h5['config']['y_scan_mv'].value

x_mv_per_deg = 1900.
y_mv_per_deg = 1900.

x_scan_deg = x_scan_mv/x_mv_per_deg
y_scan_deg = y_scan_mv/y_mv_per_deg

cone_density,cone_row_spacing = cdi.get_density_and_rowspacing(nt,si)

x_scan_m = x_scan_deg*M_PER_DEG
y_scan_m = y_scan_deg*M_PER_DEG

x_pixel_size = x_scan_m/float(n_fast)
y_pixel_size = y_scan_m/float(n_slow)
z_pixel_size = Sz

projection_keys = h5['projections'].keys()
profile_labels = h5['model']['labels'].keys()
model_profile = h5['model']['profile'][:]

pixels_per_cone = cone_row_spacing/x_pixel_size



avol = np.abs(h5['flattened_data'][vidx,:,:,:])
avol = np.transpose(avol,(0,2,1))
mprofile = h5['model']['profile'][:]
avol = avol[10:,25:,:len(mprofile)]

for z in range(avol.shape[2]):
    avol[:,:,z] = utils.gaussian_convolve(avol[:,:,z],0.5)

profile = np.mean(np.mean(avol,axis=1),axis=0)
offset = utils.translation1(profile,mprofile)[0]

costidx = h5['model']['labels']['COST'].value + offset
isosidx = h5['model']['labels']['ISOS'].value + offset
rad = 2



isos = np.mean(avol[:,:,isosidx-rad:isosidx+rad+1],axis=2)
cost = np.mean(avol[:,:,costidx-rad:costidx+rad+1],axis=2)

cones = isos+cost#np.mean(avol[:,:,isosidx-rad:costidx+rad+1],axis=2)
cones = cost

clim = np.percentile(cones,(5,99))


scones = utils.gaussian_convolve(cones,sigma=0.25,mode='same')#[25:85,45:105]
#plt.imshow(scones,cmap='gray',interpolation='none',clim=clim)

def odd(num):
    assert num==round(num)
    return num%2==1

nsize = np.floor(pixels_per_cone)
if not odd(nsize):
    nsize = nsize+1
    
cx,cy=utils.find_cones(scones,nsize,nstd=2.0,do_plot=True)
plt.figure()
plt.imshow(scones,cmap='gray',interpolation='none',clim=clim)
plt.colorbar()
plt.autoscale(False)
plt.plot(cx,cy,'b.')
plt.show()

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
nstd = 0
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
fracs = [0.1,1.0,0.1]

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

newvol = newvol + avol/4.0

nisos = np.mean(newvol[:,:,isosidx-rad:isosidx+rad+1],axis=2)
ncost = np.mean(newvol[:,:,costidx-rad:costidx+rad+1],axis=2)
ncones = np.mean(newvol[:,:,isosidx-rad:costidx+rad+1],axis=2)
    

plt.imshow(ncones,cmap='gray')#,interpolation='none')
plt.autoscale(False)
#plt.plot(cx,cy,'b.')
plt.show()
    

scatter_field = mlab.pipeline.scalar_field(newvol)
vmin,vmax = np.percentile(ncones,[5,95])
mlab.pipeline.volume(scatter_field,vmin=vmin,vmax=vmax,color=(1.0,1.0,1.0))
mlab.pipeline.image_plane_widget(scatter_field,
                                 plane_orientation='z_axes',
                                 slice_index=costidx,
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
