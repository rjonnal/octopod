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


if False:
    # in this section, the peaks are manually compared with Ocean Optics specs
    # in order to pair the expected peaks in ocean_optics_argon_lines.txt (750.387...)
    # with the peaks in the acquired image.
    
    # load the known peak locations:
    argon_peaks = np.loadtxt('./ocean_optics_argon_lines.txt')*1e-9
    argon_peaks = argon_peaks[np.where(np.logical_and(argon_peaks>np.min(L),argon_peaks<np.max(L)))]

    plt.figure()
    plt.plot(L,prof)
    plt.figure()
    plt.plot(prof)
    plt.plot(argon_peaks,np.ones(len(argon_peaks))*np.max(prof),'ks')
    plt.show()

    # using these figures, identify the wavelengths of the visible peaks (using ocean_optics_argon_lines.txt),
    # and then record, in a new file ocean_optics_argon_lines_pixel_numbers.txt, the pixel numbers associated
    # with each peak.

peak_pixels = np.loadtxt('ocean_optics_argon_lines_pixel_numbers.txt',delimiter=',')

wavelength = peak_pixels[:,0]*1e-9
pixel = peak_pixels[:,1]
valid = np.where(pixel>-1)
wavelength = wavelength[valid]
pixel = pixel[valid]

p1 = np.polyfit(pixel,wavelength,1)
p2 = np.polyfit(pixel,wavelength,2)

pixels = np.arange(2048)
p1_fit = np.polyval(p1,pixels)
p2_fit = np.polyval(p2,pixels)

plt.plot(pixels,p1_fit,'k--')
plt.plot(pixels,p2_fit,'b:')
plt.plot(pixel,wavelength,'rs')

plt.show()
