import numpy as np
from matplotlib import pyplot as plt
import sys,os


if True:

    sy = 100
    sx = 100

    im = np.random.randn(sy,sx)

    ac = np.abs(np.fft.ifft2(np.fft.fft2(im)*np.conj(np.fft.fft2(im))))
    ac = np.abs(np.fft.ifft2(np.fft.fft2(im)*np.fft.fft2(im).conjugate()))
    print ac.max()


if False:
    # a simple model of finding a line in a matrix
    # we want to know the coordinates of tar in ref
    # in this version use no resampling
    yc = 20
    xc = -40
    width = 500
    height = 500
    
    ref = np.random.randn(height,width)
    tar = np.random.randn(width)

    if xc>=0:
        tar[:-xc] = ref[yc,xc:]
    else:
        tar[-xc:] = ref[yc,:xc]
        
    ref_f = np.fft.fft2(ref)
    tar_f = np.fft.fft(tar)

    ref_fc = ref_f.conjugate()
    tar_fc = tar_f.conjugate()

    ref_ac = np.abs(np.fft.ifft2(ref_f*ref_fc))
    tar_ac = np.abs(np.fft.ifft(tar_f*tar_fc))

    joint_ac = np.sqrt(ref_ac.max())*np.sqrt(tar_ac.max())
    
    centered_xc = np.fft.fftshift(np.abs(np.fft.ifft2(ref_f*tar_fc))/joint_ac,axes=1)*np.sqrt(height)
    
    cpeaky,cpeakx = np.where(centered_xc==centered_xc.max())
    cpeaky = float(cpeaky[0])
    cpeakx = float(cpeakx[0])
    peakx = cpeakx - width // 2
    peaky = cpeaky

    print peakx,peaky,xc,yc
    
if False:
    # a simple model of finding a line in a matrix
    # we want to know the coordinates of tar in ref
    # now with oversampling

    oversample_factor = 2
    
    yc = 20
    xc = 34
    width = 532
    height = 566
    
    ref = np.random.randn(height,width)
    tar = np.random.randn(width)

    if xc>0:
        tar[:-xc] = ref[yc,xc:]
    elif xc<0:
        tar[-xc:] = ref[yc,:xc]
    else:
        tar = ref[yc,:]

    tars = [tar]

    # start cross-correlation:
    sy,sx = ref.shape
    Ny = sy * oversample_factor
    Nx = sx * oversample_factor
    
    ref_f = np.fft.fft2(ref)
    ref_fc = ref_f.conjugate()
    ref_ac = np.abs(np.fft.ifft2(ref_f*ref_fc))

    for tar in tars:
        tar_f = np.fft.fft(tar)
        tar_fc = tar_f.conjugate()
        tar_ac = np.abs(np.fft.ifft(tar_f*tar_fc))

        joint_ac = np.sqrt(ref_ac.max())*np.sqrt(tar_ac.max())
    
        centered_xc = np.fft.fftshift(np.abs(np.fft.ifft2(ref_f*tar_fc,s=(Ny,Nx)))/joint_ac,axes=1)*np.sqrt(sy)*oversample_factor**2
        xc_max = centered_xc.max()
        
        cpeaky,cpeakx = np.where(centered_xc==xc_max)
        cpeaky = float(cpeaky[0])
        cpeakx = float(cpeakx[0])

        peakx = (cpeakx - Nx // 2) / oversample_factor
        peaky = cpeaky / oversample_factor
        
        print peakx,peaky,xc,yc,xc_max
    
                   

sys.exit()

x = np.linspace(-5*np.pi,5*np.pi,64)

def noise(sig,factor=1.0):
    return np.random.randn(len(sig))*np.sqrt(np.abs(sig))*factor

def sfunc(x):
    return np.sin(x)*x**2
    #return np.random.rand(len(x))

a = sfunc(x)
b = sfunc(x+2.5)

a = a + noise(a)
b = b + noise(b)

a = (a - a.mean())/a.std()
b = (b - b.mean())/b.std()

a_f = np.fft.fft(a)
a_fc = a_f.conjugate()
a_autocorr = np.abs(np.fft.ifft(a_f*a_fc))
a_acmax = np.max(a_autocorr)

b_f = np.fft.fft(b)
b_fc = b_f.conjugate()
b_autocorr = np.abs(np.fft.ifft(b_f*b_fc))
b_acmax = np.max(b_autocorr)

num = a_f*b_fc
num = np.fft.fftshift(np.abs(np.fft.ifft(np.fft.fftshift(num))))
denom = np.sqrt(a_acmax)*np.sqrt(b_acmax)

xc = num/denom

plt.plot(x,xc)
plt.show()
sys.exit()
