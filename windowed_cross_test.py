import numpy as np
from matplotlib import pyplot as plt

image = np.zeros((3,3))
image[1,1] = 1.0

image = np.array([[0,1,0],
                  [1,0,1],
                  [0,1,0]])

#image = np.random.rand(3,3)
def printx(mat,matname):
    print matname
    sy,sx = mat.shape
    for y in range(sy):
        for x in range(sx):
            print '%0.2f %0.2f\t'%(np.real(mat[y,x]),np.imag(mat[y,x])),
        print
    print


def compute_for_sigma(sigma):
    x = np.array([-1,0,1])
    g = np.exp((-x**2)/(2*sigma**2))/np.sqrt(2*sigma**2*np.pi)
    g = g/np.max(g)

    print sum(g)
    
    ref = image
    tar = image*g


    f0 = np.fft.fft2(tar)
    f1 = np.fft.fft2(ref)
    f1c = f1.conjugate()

    num = f0*f1c
    denom = np.abs(f0)*np.abs(f1c)

    frac = num/denom
    
    complex_ir = np.fft.ifft2(frac)
    abs_ir = np.abs(complex_ir)

    printx(tar,'tar')
    printx(ref,'ref')
    printx(f0,'f0')
    printx(f1c,'f1c')
    printx(num,'num')
    printx(denom,'denom')
    printx(frac,'frac')
    printx(abs_ir,'abs_ir')

compute_for_sigma(500.0)
compute_for_sigma(0.5)
