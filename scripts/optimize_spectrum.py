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


# previous values
L0_original, dL_original = (9.2633165829145736e-07, -9.4070351758793966e-11)

# realistic bounds:
L0_limits = (930e-9,1060e-9)
dL_limits = (-10e-11,-11e-11)

# crazy bounds
L0_limits = (530e-9,1560e-9)
dL_limits = (-100e-11,-5e-11)

# liberal bounds
L0_limits = (860e-9,1100e-9)
dL_limits = (-15e-11,-8e-11)

# bounds based on previous values
L0_limits = (L0_original*.7,L0_original*1.3)
dL_limits = (dL_original*1.3,dL_original*.7)

Ns = (20,20)

output = np.zeros(Ns)

range1 = np.linspace(L0_limits[0],L0_limits[1],Ns[0])
range2 = np.linspace(dL_limits[0],dL_limits[1],Ns[1])

max_val = -np.inf
for idx1,L0 in enumerate(range1):
    print idx1
    for idx2,dL in enumerate(range2):
        val = obj(L0,dL)
        if val>max_val:
            winners = L0,dL
            max_val = val
        output[idx1,idx2] = val

print winners
plt.figure()
plt.imshow(output,interpolation='none',aspect='auto')
plt.autoscale(False)
ymax,xmax = np.where(output==np.max(output))
ph = plt.plot(xmax,ymax,'ks')[0]
ax = plt.gca()
#ax.set_xticklabels(range2[::5])
#ax.set_yticklabels(range1[::5])
plt.colorbar()

def save_trial():
    L00 = '%0.1f'%(L0_limits[0]*1e9)
    L01 = '%0.1f'%(L0_limits[1]*1e9)
    dL0 = '%0.1f'%(dL_limits[0]*1e11)
    dL1 = '%0.1f'%(dL_limits[1]*1e11)
    L0N = '%d'%(Ns[0])
    dLN = '%d'%(Ns[1])
    w0 = '%0.1f'%(winners[0]*1e9)
    w1 = '%0.1f'%(winners[1]*1e11)
    outfn = './optimize_spectrum_output/%s_%s_%s_%s_%s_%s_%s_%s.png'%(L00,L01,dL0,dL1,L0N,dLN,w0,w1)
    plt.savefig(outfn,dpi=300)

save_trial()

plt.show()
