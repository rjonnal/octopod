# Paths to raw data; all valid paths will be used by the
# file manager.  For each key/value pair in the data_paths
# dictionary, the key will be used to label data sets in
# the processing steps.
data_paths = {}
#data_paths['2g_aooct'] = 'D:/rjonnal/Dropbox/Share/2g_aooct_data/Data/'
data_paths['2g_aooct'] = '/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data'
data_paths['hroct'] = '/home/rjonnal/data/Dropbox/Share/hroct_data/Data'

db_file = '/home/rjonnal/db/octopod.db'

raw_data_extension = 'unp'
raw_bytes_per_pixel = 2

source_spectra = {}
source_spectra['hroct'] = {'L0':801e-9,'dL':6.2e-11}
source_spectra['2g_aooct'] = {'L0':920e-9,'dL':-9.5e-11}

dc_cutoff = 20
