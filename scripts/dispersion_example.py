from octopod import *
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)

fm = FileManager()
files = fm.get('hroct')

for fn in files:
    print fn
    d = Dataset(fn)
    h5 = d.get_h5_handle()

    do = DispersionOptimizer(h5)
    #tf = do.make_test_frame(2000)
    tf = h5['raw_data'][0][10:25:5]
    tf = np.reshape(tf,(3000,2048))
    dccoef = do.optimize(tf)
    h5.close()

