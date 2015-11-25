NUM_PROCESSORS = 1

# Server information
SERVER_USERNAME = 'vsri'
SERVER_PASSWORD = 'Augen'
SERVER_IP = '152.79.39.148'
SERVER_WWW_ROOT = '/home/vsri/www/data/'


# DATA_ROOTS is a list of one or more paths were OCT data is to be
# found on the acquisition machine. Each item in it should be an
# absolute path to a folder containing experiment date folders
# (e.g. 2012.01.05 or 2_10_2013).

DATA_ROOTS = ['D:/afocal_AO-OCT/Data/']
# DATA_ROOTS = ['/dose/Share/afocal_aooct_data/Data/','/dose/Share/afocal_aooct_data/Data_Unused/','/home/rjonnal/Dropbox/Share/afocal_aooct_data/Data/']
#DATA_ROOTS = ['H:/AOOCT_Data/2014/']

# MASTER LOG
MASTER_LOG = 'D:/afocal_AO-OCT/Data/masterlog.txt'
# MASTER_LOG = '/home/rjonnal/data/Dropbox/Share/afocal_aooct_data/Data/masterlog.txt'

# specify the end points of the source spectrum:

# MAPPING_POLYNOMIAL contains coefficients for an Nth order polynomial
# L(x) = MP[0]*x^N + MP[1]*x^(N-1) + ... + MP[-1], such that L(x) describes
# the wavelength (in meters) of the light falling on the xth pixel of the sensor.
MAPPING_POLYNOMIAL = [6.751343e-11,7.747e-7] # linear polynomial for 2g-ao-oct system, derived from ocean optics values
#MAPPING_POLYNOMIAL = [ -2.67600024e-15, 7.30246168e-11, 7.83405956e-07] # quadratic polynomial for 1g-ao-oct system, derived from ocean optics values
# MAPPING_POLYNOMIAL = [0.065e-9, 807e-09] # bioptigen values

# is the image flipped (i.e. are we doing EDI?):
FLIPPED = True
FAST_MIN = 0

DC_REMOVAL_LENGTH = 20
DEPTH_MIN = 1024 + DC_REMOVAL_LENGTH
DEPTH_MAX = 1500

#DISPERSION_COEFS = [-6e-18,-2.2e-11,0.,0.] # 2G-AO-OCT working values
DISPERSION_COEFS = [0.0,0.0,0.,0.] # 1G-AO-OCT working values
DISPERSION_COEFS_UBOUNDS = [4e-16,-2e-10,0.0,0.0] # 1G-AO-OCT working values
DISPERSION_COEFS_LBOUNDS = [2.5e-16,-4e-10,0.0,0.0] # 1G-AO-OCT working values


# technical parameters:
# the number of pixels, symmetrically arranged about the DC peak, to zero in order to eliminate the DC
# theoretically this should be the PSF, but in practice it can be much larger due to dispersion
NOISE_FLOOR_MAX_SLOPE = 30


# the number of standard deviations (of the noise) above the noise floor that should indicate
# the presence of signal
# if the noise were normally distributed, 4 would probably suffice for a 2k vector, but 5 or 6
# might be safer
SIGNAL_THRESHOLD_N_STD = 5

# when automatically cropping out the retina, what's the maximum full retinal thickness to crop
RETINAL_THICKNESS_M = 500e-6

# when using model-based localization of a-line peaks, how far from the model location should we look?
PEAK_SEARCH_RADIUS = 4

N_EYE = 1.38
N_OS = 1.41
N_AIR = 1.0
N_GLASS = 1.5
N_LIPID = 1.48
N_WATER = 1.33


##########################################################################################
# OBSOLETE / DEPRECATED SETTINGS:
# standard small volume settings:
#DEFAULT_NV = 3
#DEFAULT_NX = 150
#DEFAULT_NY = 170
#DEFAULT_NZ = 2048

# alternate settings for larger volume:
#DEFAULT_NV = 1
#DEFAULT_NX = 400
#DEFAULT_NY = 500
#DEFAULT_NZ = 2048


