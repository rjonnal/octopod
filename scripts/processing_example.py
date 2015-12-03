from octopod import *
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)

fm = FileManager()
files = fm.get('hroct')

for fn in files:
    d = Dataset(fn)
    h5 = d.get_h5_handle()
    p = OCTProcessor(h5)
    p.cleanup()
    p.run()
    h5.close()


