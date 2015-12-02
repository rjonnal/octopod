from octopod import *
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)

fm = FileManager()
files = fm.get('hroct')

acf = AcquisitionParameterFile()
for fn in files:
    acf.make_xml_file(fn)


