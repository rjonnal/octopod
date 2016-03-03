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

plt.figure()
plt.plot(L,np.mean(im,axis=0))
plt.show()

