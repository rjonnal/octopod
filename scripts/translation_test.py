import numpy as np
from matplotlib import pyplot as plt

target = np.zeros((10,10))
template = np.zeros((10,10))

target[1,:] = 1.0
template[0,:] = 1.0

ftar = np.fft.fft2(target)
ftem = np.fft.fft2(template)

xc = np.abs(np.fft.ifft2((ftar * np.conj(ftem)) / (np.abs(ftar) * np.abs(ftem))))

plt.subplot(3,2,1)
plt.imshow(target,interpolation='none',aspect='auto')
plt.subplot(3,2,2)
plt.imshow(template,interpolation='none',aspect='auto')
plt.subplot(3,2,3)
plt.imshow(np.abs(ftar),interpolation='none',aspect='auto')
plt.subplot(3,2,4)
plt.imshow(np.abs(ftem),interpolation='none',aspect='auto')
plt.subplot(3,2,5)
plt.imshow(xc,interpolation='none',aspect='auto')

plt.show()
