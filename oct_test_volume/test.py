from octopod import *
import sys

ds = Dataset('./oct_test_volume_2T.unp')
ds.initialize('2g_aooct')
ds.optimize_dispersion()
sys.exit()
ds.process()

frame = ds.h5.get('processed_data')[0,0,:,:]
plt.imshow(np.abs(frame),interpolation='none',aspect='auto')
plt.show()
