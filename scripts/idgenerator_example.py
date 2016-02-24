from octopod.H5Utils import *
import glob,os
import sys

idg = IDGenerator()

eg = EccentricityGuesser()

subjects = ['Jonnal_Ravi','Werner_Max']
roots = ['/home/rjonnal/Dropbox/Share/2g_aooct_data/Data/2016.01.29','/home/rjonnal/Dropbox/Share/2g_aooct_data/Data/2016.02.02']
dates = ['20160129','20160202']

for subject,root,date in zip(subjects,roots,dates):
    flist = glob.glob(os.path.join(root,'*.hdf5'))
    for filename in flist:
        idg.create_ids(filename,subject,date)
        eg.guess_to_h5(filename)
