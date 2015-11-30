from octopod import *
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)

fm = FileManager()
fn = fm.get_one('hroct',['carmen','line_3'])
ds = Dataset(fn)

test = {}
test['time']={'Data_Acquired_at':'11/17/2015 4:21:01 PM'}
test['volume_size']={'width':'2048','height':'400','number_of_frames':'500','number_of_volumes':'1'}
test['scanning_parameters']={'x_scan_range':'3168','x_scan_offset':'0','y_scan_range':'2112','y_scan_offset':'0','number_of_bm_scans':'1'}
test['dispersion_parameters']={'c2':'1.35e-5','c3':'0'}
