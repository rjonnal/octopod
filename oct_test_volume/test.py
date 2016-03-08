from octopod import *

ds = Dataset('./oct_test_volume_2T.unp')
ds.initialize('2g_aooct')
ds.optimize_dispersion()
ds.process()
