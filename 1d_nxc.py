import numpy as np

a = np.zeros((10))
b = np.zeros((10))

a[5] = 1
b[6] = 1

af = np.fft.fft(a)
bf = np.fft.fft(b)
