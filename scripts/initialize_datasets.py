from octopod import *
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)

fm = FileManager()
files = fm.get('hroct')

for fn in files:
    d = Dataset(fn)
    d.initialize()
    sys.exit()
