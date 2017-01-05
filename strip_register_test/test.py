from octopod import utils
import numpy as np
from scipy.misc import imread

im = imread('kids.jpg')

im1 = im[5:,5:]
im2 = im[:-5,:-5]


yp,xp,g = utils.strip_register(im1,im2,1,1.0,True)

