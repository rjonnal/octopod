import numpy as np
from matplotlib import pyplot as plt
import sys,os



if True:
    # a simple model of finding a line in a matrix
    # we want to know the coordinates of tar in ref
    
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
    plt.imshow(centered_xc,cmap='gray',interpolation='none')
    plt.colorbar()
    plt.show()
    
    
    cpeaky,cpeakx = np.where(centered_xc==centered_xc.max())
    cpeaky = float(cpeaky[0])
    cpeakx = float(cpeakx[0])
    peakx = cpeakx - width // 2
    peaky = cpeaky

    print peakx,peaky
    
                   

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
