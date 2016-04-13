import numpy as np
from matplotlib import pyplot as plt

a = np.zeros((10))
b = np.zeros((10))

a[5] = 1
b[6] = 1

af = np.fft.fft(a)
bfc = np.conj(np.fft.fft(b))

prod = af*bfc
plt.plot(np.imag(prod))
plt.show()


