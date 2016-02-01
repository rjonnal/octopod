from octopod import *
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)
import sys,os
import shutil

system_label = '2g_aooct'
fm = FileManager()

files = fm.get(system_label,['2016.01.29'])

home_dir = 'D:/rjonnal/Dropbox/Share/2g_aooct_data/Data/2016.01.29'

for fn in files:
    d = Dataset(fn)
    h5fn = d.h5fn
    home_fn = os.path.join(home_dir,os.path.split(h5fn)[1])
    print home_dir
    # try:
        # del h5['dispersion/coefficients']
    # except:
        # pass
    # h5.require_group('dispersion')
    # h5['dispersion'].create_dataset('coefficients',data=np.array([0.0,0.0,0.0,0.0]))
    # do = DispersionOptimizer(h5)
    # tf = do.make_test_frame(1000)
    # dccoef = do.optimize(tf)
    # p = OCTProcessor(h5)
    # p.cleanup()
    # p.run()
    # h5.close()
    shutil.copyfile(h5fn,home_fn)
    


