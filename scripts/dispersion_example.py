from octopod import *
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)

fm = FileManager()
fn = fm.get_one('hroct',['carmen','line_1'])

d = Dataset(fn)
h5 = d.get_h5_handle()

dccoefs = [0.0,0.0,0.0,0.0]
do = DispersionOptimizer(h5)
inf = h5['raw_data'][0][0]
outf1 = do.process_frame(inf,[0.0,0.0,0.0,0.0])
outf2 = do.process_frame(inf,[-4.5999999999999996e-17, -1.1999999999999991e-11, 0.0, 0.0])

plt.figure()
plt.imshow(np.abs(outf1[700:900,:200]))
plt.colorbar()
plt.figure()
plt.imshow(np.abs(outf2[700:900,:200]))
plt.colorbar()
plt.show()


sys.exit()
tf = do.make_test_frame(50)
dccoef = do.get(tf)
h5.close()

print dccoef
