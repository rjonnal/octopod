from octopod import *
import numpy as np

h5 = H5('test.hdf5')

arr = np.array([1,2,3])

h5.put('arr',arr,short_descriptor='a vector')

val = 10.0

h5.put('blah/val',val,short_descriptor='a scalar')

def info(item):
    print type(item)
    print item
    
info(h5.get('/arr'))
info(h5.get('blah/val'))




