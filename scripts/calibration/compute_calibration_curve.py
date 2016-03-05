import numpy as np
import scipy as sp
from PIL import Image
from matplotlib import pyplot as plt
import sys
from octopod import octopod_config as ocfg

fn_root = '2g_aooct_oceanoptics_calibration'
ext = 'tif'

imfn = '%s.%s'%(fn_root,ext)

im = np.array(Image.open(imfn))

# ballpark lambda
L = np.arange(2048)*ocfg.source_spectra['2g_aooct']['dL'] + ocfg.source_spectra['2g_aooct']['L0']


im = im - np.min(im)

thresh = 140
prof = np.mean(im,axis=0)

#prof[np.where(prof<thresh)]=0


# load the known peak locations:
argon_peaks = np.loadtxt('./ocean_optics_argon_lines.txt')*1e-9
argon_peaks = argon_peaks[np.where(np.logical_and(argon_peaks>np.min(L),argon_peaks<np.max(L)))]

plt.figure()
plt.plot(L,prof)
plt.plot(argon_peaks,np.ones(len(argon_peaks))*np.max(prof),'ks')
    


plt.show()

