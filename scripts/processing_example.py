from octopod import *
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)

system_label = '2g_aooct'

fm = FileManager()
files = fm.get(system_label,['Edric'])

for fn in files:
    d = Dataset(fn)
    d.initialize(system_label)
    h5 = d.get_h5_handle()
    try:
        del h5['dispersion/coefficients']
    except:
        pass
    h5.require_group('dispersion')
    h5['dispersion'].create_dataset('coefficients',data=np.array([0.0,0.0,0.0,0.0]))
    #do = DispersionOptimizer(h5)
    #tf = do.make_test_frame(1000)
    #dccoef = do.optimize(tf)
    p = OCTProcessor(h5)
    p.cleanup()
    p.run()
    h5.close()



