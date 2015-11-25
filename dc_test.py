import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import hilbert
import sys


x = np.arange(100)
#x = x - x.mean()

dcrand = np.random.rand(10)
vec = np.ones((len(x),10))*dcrand

vec = hilbert(vec,axis=0)
d = np.exp(1j*np.polyval([1.0e-7,1.0e-3,0.0,0.0],x))
df = np.fft.fftshift(np.abs(np.fft.fft(d)))
vec = (vec.T * d).T
vec = np.fft.fftshift(np.fft.fft(vec,axis=0))


plt.imshow((np.abs(vec).T/df).T)
plt.show()


