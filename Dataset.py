import sys,os,time
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import logging
import octopod_config as ocfg

from octopod.AcquisitionParameterFile import AcquisitionParameterFile
from octopod.DataStore import H5
from octopod.Misc import EccentricityGuesser, IDGenerator
from octopod.Processor import OCTProcessor
from octopod.DispersionOptimizer import DispersionOptimizer
from octopod.Cropper import Cropper
from octopod.Model import Model
from octopod.BScanAligner import BScanAligner
from octopod.Reporter import Reporter

class Dataset:

    def __init__(self,raw_data_filename):
        """Initiate a Dataset object based on a .unp raw OCT data file."""
        self.raw_data_filename = raw_data_filename
        self.logger = logging.getLogger(__name__)
        self.logger.info('Creating Dataset based on %s.'%raw_data_filename)
        raw_temp,raw_extension = os.path.splitext(self.raw_data_filename)
        self.h5fn = raw_temp + '.hdf5'
        self.xml_fn = raw_temp + '.xml'
        self.h5 = H5(self.h5fn)

    def get_h5_handle(self):
        return self.h5

    def ecc_to_h5(self):
        eg = EccentricityGuesser()
        si_ecc,nt_ecc = eg.guess(self.raw_data_filename)
        self.h5.delete('eccentricity')
        self.h5.put('eccentricity/superior_inferior',si_ecc)
        self.h5.put('eccentricity/nasal_temporal',nt_ecc)
        self.h5.put('eccentricity/superior_and_nasal_are_negative',[np.nan])
        

    def get_acquisition_date(self):
        time_struct = time.localtime(os.path.getmtime(self.raw_data_filename))
        return '%04d%02d%02d'%(time_struct.tm_year,time_struct.tm_mon,time_struct.tm_mday)
        
    def make_ids(self):
        idg = IDGenerator()
        at = self.get_acquisition_date()
        try:
            head,tail = os.path.split(self.raw_data_filename)
            sfn = os.path.join(head,'subject_name.txt')
            fid = open(sfn)
            sn = fid.readline().strip()
            fid.close()
            sn = sn.replace(', ','_')
            sn = sn.replace(',','_')
        except Exception as e:
            self.logger.warning('No subject_name.txt file found. Using doe_john.')
            sn = 'doe_john'
            
        idg.create_ids(self.h5,sn,at)

    def process(self):
        op = OCTProcessor(self.h5)
        op.run()

    def optimize_dispersion(self,n_lines=1000):
        do = DispersionOptimizer(self.h5)
        do.optimize(n_lines=n_lines)


    def flip(self):
        test = self.h5['processed_data'][0,0,:,:]
        test = np.abs(test)
        test = np.mean(test,axis=1)
        plt.plot(test)
        plt.show()
        sys.exit()

    def crop(self):
        c = Cropper(self.h5)
        c.crop()

    def model(self):
        m = Model(self.h5)
        m.click_label()

    def align(self,do_plot=False):
        bsa = BScanAligner(self.h5)
        bsa.align_volumes(do_plot=do_plot)


    def label(self):
        m = Model(self.h5)
        m.write_axial_alignment()
        m.write_volume_labels()


    def report(self,do_log=False):
        report_directory = self.h5.filename.replace('.hdf5','')+'_report'
        r = Reporter(self.h5,report_directory=report_directory)
        try:
            r.processed_report(do_log=do_log)
        except Exception as e:
            self.logger.info('processed_report: Could not complete: %s'%e)
        try:
            r.model_report()
        except Exception as e:
            self.logger.info('model_report: Could not complete: %s'%e)
        try:
            r.dispersion_report()
        except Exception as e:
            self.logger.info('dispersion_report: Could not complete: %s'%e)
        try:
            r.projections_report()
        except Exception as e:
            self.logger.info('projections_report: Could not complete: %s'%e)


            
    def initialize(self,system_label):

        # write parameters from the XML file to the h5 file
        apf = AcquisitionParameterFile()
        apf.translate_xml_to_h5(self.xml_fn,self.h5)
        
        n_vol = self.h5.get('/config/n_vol')[()]
        n_slow = self.h5.get('/config/n_slow')[()]
        n_fast = self.h5.get('/config/n_fast')[()]
        n_depth = self.h5.get('/config/n_depth')[()]

        raw_store = self.h5.make('raw_data',(n_vol,n_slow,n_fast,n_depth),dtype='u2')

        with open(self.raw_data_filename,'rb') as fid:
            for vol_index in range(n_vol):
                position = vol_index * n_depth * n_fast * n_slow * ocfg.raw_bytes_per_pixel
                fid.seek(position,0)
                vol = np.fromfile(fid,dtype=np.uint16,count=n_slow*n_depth*n_fast)
                vol = vol.reshape(n_slow,n_fast,n_depth)
                raw_store[vol_index,:,:,:] = vol

        p = ocfg.source_spectra[system_label] # polynomial for generating dataset's pixel->lambda function
        
        self.L = np.polyval(p,np.arange(n_depth))
        self.k_in = (2.0*np.pi)/self.L
        self.k_out = np.linspace(self.k_in[0],self.k_in[-1],n_depth)

        self.h5.put('L',self.L)
        self.h5.put('k_in',self.k_in)
        self.h5.put('k_out',self.k_out)
        self.h5.put('system_label',system_label)

        #self.h5.put('dispersion/coefficients',[0.0,0.0])
        
        self.ecc_to_h5()
        self.make_ids()
                
def test():
    ds = Dataset('./oct_test_volume/oct_test_volume_2T.unp')
    ds.initialize('2g_aooct')
        
if __name__=='__main__':
    test()
