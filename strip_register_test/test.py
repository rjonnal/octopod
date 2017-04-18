from octopod import utils
import numpy as np
from scipy.misc import imread
from matplotlib import pyplot as plt

#im = np.load('../../raster_tools/images/grass.npy')
im = imread('./kids.jpg')

ref = im[200:500,100:400]
tar = im[201:501,102:402]

plt.figure()
plt.subplot(1,2,1)
plt.imshow(ref,interpolation='none',cmap='gray')
plt.title('reference')
plt.subplot(1,2,2)
plt.imshow(tar,interpolation='none',cmap='gray')
plt.title('target')
plt.show()

yp,xp,g = utils.strip_register(tar,ref,5,10.0,False)

plt.figure()
plt.subplot(2,1,1)
plt.plot(xp)
plt.plot(yp)
plt.subplot(2,1,1)
plt.plot(g)
plt.show()

# output should have plots of y=-1 and x=-2, revealing amounts
# by which TAR should be shifted to match REF
