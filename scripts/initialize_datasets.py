from octopod import *
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)

fm = FileManager()
system_label = 'hroct'
files = fm.get(system_label)

for fn in files:
    d = Dataset(fn)
    d.initialize(system_label)
