from matplotlib import pyplot as plt
import numpy as np
from octopod import *
import glob

unp_list = glob.glob('*.unp')

for unp in unp_list:
    #ds = Dataset(unp)
    #ds.initialize('2g_aooct')
    #ds.optimize_dispersion(n_lines=1000)

    h5fn = unp.replace('.unp','')+'.hdf5'
    h5 = H5(h5fn,mode='w')

    if False:
        p = OCTProcessor(h5)
        p.run()


    if False:
        c = Cropper(h5)
        c.crop()

    if False:
        pd = h5['processed_data']
        rd = h5['raw_data']
        print (np.prod(pd.shape)*8+np.prod(rd.shape)*2)/1024./1024.

    if False:
        m = Model(h5)
        m.click_label()

    if True:
        m = Model(h5)
        m.write_axial_alignment()
        m.write_volume_labels()
    
    #r = Reporter(h5,report_directory='./%s'%os.path.splitext(unp)[0])
    
    #r.dispersion_report()

    
