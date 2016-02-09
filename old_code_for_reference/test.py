from octopod import *
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)

fm = FileManager()
fn = fm.get_one('hroct',['carmen','line_3'])
ds = Dataset(fn)

acf = AcquisitionParameterFile()
acf.make_xml_file(fn)