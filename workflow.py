from octopod import *
import glob

unp_fn_list = glob.glob('./*.unp')

stage = 1

if stage==1:
    # one loop for non-interactive, fully automated steps:
    for unp_fn in unp_fn_list:
        # make a dataset and initialize it, specifying the system with which the data was collected
        ds = Dataset(unp_fn)
        ds.initialize('2g_aooct')
        #ds.initialize('hroct')

        # call the dispersion optimizer
        # this is equivalent to instantiating a DispersionOptimizer object and running
        # its optimizer; instead we use the shortcut in the Dataset object
        # the default value of n_lines is 1000; here we use only 100 to speed things
        # up, for demonstration
        ds.optimize_dispersion(n_lines=100)


        # call the processor
        # this is equivalent to instantiating a Processor object and running it
        ds.process()

if stage==2:
    # a second loop for interactive steps:
    for unp_fn in unp_fn_list:
        ds = Dataset(unp_fn)
        ds.initialize('2g_aooct')

        # interactive step:
        # now that we have a complete processed data cube, let's crop it to
        # save on storage and ram; again, this is a shortcut to the Cropper class
        # and its crop method:
        ds.crop()

        # make and store aligned versions of the volumes:
        ds.align(do_plot=True)

        # interactive step:
        # make a model and label its features:
        ds.model()

if stage==3:
    # a third loop for automated steps:
    for unp_fn in unp_fn_list:
        ds = Dataset(unp_fn)
        ds.initialize('2g_aooct')
        ds.label()

        ds.report()
