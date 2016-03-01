import sys,os
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import logging
from octopod.AcquisitionParameterFile import AcquisitionParameterFile
from octopod.DataStore import H5
from octopod.H5Utils import EccentricityGuesser
import octopod_config as ocfg


class Dataset:

    def __init__(self,raw_data_filename):
        """Initiate a Dataset object based on a .unp raw OCT data file."""
        self.raw_data_filename = raw_data_filename
        self.logger = logging.getLogger(__name__)
        self.logger.info('Creating Dataset based on %s.'%raw_data_filename)
        raw_temp,raw_extension = os.path.splitext(self.raw_data_filename)
        self.h5fn = raw_temp + '.hdf5'
        self.xml_fn = raw_temp + '.xml'

    def get_h5_handle(self):
        return self.h5

    def ecc_to_h5(self):
        eg = EccentricityGuesser()
        si_ecc,nt_ecc = eg.guess(self.raw_data_filename)

        self.h5.delete('eccentricity')
        self.h5.put('eccentricity/superior_infierior',si_ecc)
        self.h5.put('eccentricity/nasal_temporal',si_ecc)
        self.h5.put('eccentricity/superior_and_nasal_are_negative',[np.nan])
        

    def initialize(self,system_label):

        self.h5 = H5(self.h5fn)

        # write parameters from the XML file to the h5 file
        apf = AcquisitionParameterFile()
        apf.translate_xml_to_h5(self.xml_fn,self.h5)
        
        n_vol = self.h5['/config/n_vol'][()]
        n_slow = self.h5['/config/n_slow'][()]
        n_fast = self.h5['/config/n_fast'][()]
        n_depth = self.h5['/config/n_depth'][()]

        raw_store = self.h5.make('raw_data',(n_vol,n_slow,n_fast,n_depth),dtype='u2')

        with open(self.raw_data_filename,'rb') as fid:
            for vol_index in range(n_vol):
                position = vol_index * n_depth * n_fast * n_slow * ocfg.raw_bytes_per_pixel
                fid.seek(position,0)
                vol = np.fromfile(fid,dtype=np.uint16,count=n_slow*n_depth*n_fast)
                vol = vol.reshape(n_slow,n_fast,n_depth)
                raw_store[vol_index,:,:,:] = vol

        L0 = ocfg.source_spectra[system_label]['L0']
        dL = ocfg.source_spectra[system_label]['dL']
        self.L = np.arange(n_depth)*dL+L0
        self.k_in = (2.0*np.pi)/self.L
        self.k_out = np.linspace(self.k_in[0],self.k_in[-1],n_depth)

        self.h5.put('L',self.L)
        self.h5.put('k_in',self.k_in)
        self.h5.put('k_out',self.k_out)
        self.h5.put('system_label',system_label)

        self.ecc_to_h5()
        self.h5.close()
                

if __name__=='__main__':

    ds = Dataset('./oct_test_volume/oct_test_volume.unp')
    ds.initialize('2g_aooct')
