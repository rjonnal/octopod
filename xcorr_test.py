import numpy as np
from matplotlib import pyplot as plt
base = np.random.rand(100,120)


def xcorr(x,y,oversample_factor):
    xsy,xsx = x.shape
    ysy,ysx = y.shape
    assert xsy==ysy and xsx==ysx
    Ny,Nx = ysy*oversample_factor,ysx*oversample_factor
    return np.abs(np.fft.ifft2(np.fft.fft2(x)*np.fft.fft2(y).conjugate(),s=(Ny,Nx)))

def nxcorr(x,y,oversample_factor=1):
    xac = xcorr(x,x,oversample_factor)
    yac = xcorr(y,y,oversample_factor)
    xac_max = np.max(xac)
    yac_max = np.max(yac)

    xc = xcorr(x,y,oversample_factor)
    xc_max = np.max(xc)

    #print xc_max,xac_max,yac_max
    
    corr_coef = xc_max/np.sqrt(yac_max)/np.sqrt(xac_max)
    
    return corr_coef


mine = []
theirs = []
noises = []

for noise_std in np.linspace(0.0,5.0,100):
    noises.append(noise_std)
    a = base+np.random.randn(100,120)*noise_std
    a[50:,:] = 0.0
    b = base+np.random.randn(100,120)*noise_std

    theirs.append(np.corrcoef(a.ravel(),b.ravel())[0,1])
    mine.append(nxcorr(a,b,2))


mine = np.array(mine)
theirs = np.array(theirs)
noises = np.array(noises)

mine = mine**2
plt.figure()
plt.plot(noises,theirs)
plt.plot(noises,mine)
plt.figure()
plt.plot(mine,theirs)
plt.show()
