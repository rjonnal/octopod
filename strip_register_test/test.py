from octopod import utils
import numpy as np
from scipy.misc import imread
from matplotlib import pyplot as plt

im = np.load('../../raster_tools/images/grass.npy')

im1 = im[10:100,10:100]
im2 = im[11:101,12:102]

plt.figure()
plt.subplot(1,2,1)
plt.imshow(im1,interpolation='none',cmap='gray')
plt.subplot(1,2,2)
plt.imshow(im2,interpolation='none',cmap='gray')
plt.show()

yp,xp,g = utils.strip_register(im1,im2,5,10.0,True)

plt.figure()
plt.subplot(2,1,1)
plt.plot(xp)
plt.plot(yp)
plt.subplot(2,1,1)
plt.plot(g)
plt.show()
