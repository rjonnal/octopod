from octopod.Processor import process
import h5py
import numpy as np
from scipy import optimize
from matplotlib import pyplot as plt

fn = '/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.01.29/3T_1_focus30-10reps.hdf5'

h5 = h5py.File(fn)

frame = h5['raw_data'][0,10,2:,:]

def obj(L0,dL):
    L = np.arange(2048)*dL+L0
    k_in = (2.0*np.pi)/L
    k_out = np.linspace(k_in[0],k_in[-1],2048)
    proc = process(frame,k_in,k_out,[0.0,0.0,0.0,0.0])
    #return np.max(np.abs(proc))
    return np.mean(np.max(np.abs(proc),axis=0))

# realistic bounds:
L0_limits = (930e-9,1060e-9)
dL_limits = (-10e-11,-11e-11)

# crazy bounds
L0_limits = (530e-9,1560e-9)
dL_limits = (-100e-11,-5e-11)


Ns = (10,10)

output = np.zeros(Ns)

range1 = np.linspace(L0_limits[0],L0_limits[1],Ns[0])
range2 = np.linspace(dL_limits[0],dL_limits[1],Ns[1])

for idx1,L0 in enumerate(range1):
    print idx1
    for idx2,dL in enumerate(range2):
        output[idx1,idx2] = obj(L0,dL)

        
print output
plt.figure()
plt.imshow(output**2.0,interpolation='none',aspect='auto')
plt.show()
