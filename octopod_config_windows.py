# Paths to raw data; all valid paths will be used by the
# file manager.  For each key/value pair in the data_paths
# dictionary, the key will be used to label data sets in
# the processing steps.

data_root = 'D:/Data/'

data_paths = {}
data_paths['2g_aooct'] = 'D:/Data/'

db_file = 'D:/octopod.db'

raw_data_extension = 'unp'
h5_extension = 'hdf5'
raw_bytes_per_pixel = 2

source_spectra = {}

# old approach, with explicit L0 and dL terms
# source_spectra['hroct'] = {'L0':801e-9,'dL':6.2e-11}
# original values, from note on table:
# source_spectra['2g_aooct'] = {'L0':920e-9,'dL':-9.5e-11}
# values using ocean optics calibration procedure
# (see octopod/scripts/calibration for details)
# source_spectra['2g_aooct'] = {'L0':890e-9,'dL':-7.03e-11}

# new approach: assign a list of length k, [c_{k-1}, c_{k-2}, ... c_0],
# specifying the polynomial lambda = c_{k-1} x^{k-1} + c_{k-2} x^{k-2} + ... c_0,
# where the polynomial maps pixel number (x) onto wavelength (lambda).
# See numpy.polyval for identical interpretation of coefficient lists.
source_spectra['hroct'] = [6.2e-11,801e-9]
source_spectra['2g_aooct'] = [-6.79e-11,8.86e-7]
source_spectra['2g_aooct'] = [-1.82e-15,-6.31e-11,8.83e-7]


dc_cutoff = 50
model_database = 'D:/Data_Share/Dropbox/Share/global_data/axial_model/database.hdf5'
dispersion_database = 'D:/Data_Share/Dropbox/Share/global_data/dispersion/database.hdf5'

dispersion_3_max = 10.0
dispersion_3_min = -10.0
dispersion_3_multiplier = 1e-17
dispersion_3_step_size = 1e-2
dispersion_2_max = 10.0
dispersion_2_min = -10.0
dispersion_2_multiplier = 1e-11
dispersion_2_step_size = 1e-2

x_mv_per_deg = 3168.
y_mv_per_deg = 2112.
