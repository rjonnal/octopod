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
do.process_frame(dccoefs)
h5.close()


