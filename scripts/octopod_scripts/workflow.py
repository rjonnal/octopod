from octopod import *
import glob
import sys,os

dlist = glob.glob('D:/Data/2016.04.25*')

hdf5_fn_list = []
for d in dlist:
    hdf5_fn_list = hdf5_fn_list + glob.glob(os.path.join(d,'*.hdf5'))
    
lengths = []
for hdf5 in hdf5_fn_list:
    print hdf5,H5(hdf5).get_shape('processed_data')

sys.exit()
unp_fn_list = []
for d in dlist:
    unp_fn_list = unp_fn_list + glob.glob(os.path.join(d,'*.unp'))

stage = 2

# one loop for non-interactive, fully automated steps:
if stage==1:
    for unp_fn in unp_fn_list:
        # make a dataset and initialize it, specifying the system with which the data was collected
        ds = Dataset(unp_fn)
        ds.initialize('2g_aooct')

        # call the dispersion optimizer
        # this is equivalent to instantiating a DispersionOptimizer object and running
        # its optimizer; instead we use the shortcut in the Dataset object
        # the default value of n_lines is 1000; here we use only 100 to speed things
        # up, for demonstration
        ds.optimize_dispersion(n_lines=100)

        # call the processor
        # this is equivalent to instantiating a Processor object and running it
        ds.process()

# a second loop for interactive steps:
if stage==2:    
    for unp_fn in unp_fn_list:
        ds = Dataset(unp_fn)
        # interactive step:
        # now that we have a complete processed data cube, let's crop it to
        # save on storage and ram; again, this is a shortcut to the Cropper class
        # and its crop method:
        if ds.h5.get_shape('processed_data')[2]==1024:
            ds.crop()

        # make and store aligned versions of the volumes:
        #ds.align(do_plot=True)
        
        # interactive step:
        # make a model and label its features:
        #ds.model()

# a third loop for automated steps:
if stage==3:
    for unp_fn in unp_fn_list:
        ds = Dataset(unp_fn)
        ds.label()
        ds.report()
