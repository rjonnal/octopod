from octopod import Hive
import numpy as np

location = '/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2017.09.18/test_set'

h = Hive(location)
data = np.random.rand(100,100)
h.put('/grum',data,overwrite=True)
h.put('/grum/tum/mum',data)
print h.has('grum/tum/mum')
print h.get_shape('grum/tum/mum')
h.delete('grum/tum/mum')
print h.has('grum/tum/mum')




