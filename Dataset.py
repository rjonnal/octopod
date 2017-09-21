import sys,os,time
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import logging
import octopod_config as ocfg

from octopod.AcquisitionParameterFile import AcquisitionParameterFile
from octopod.DataStore import Hive
from octopod.Misc import EccentricityGuesser, IDGenerator
from octopod.Processor import OCTProcessor
from octopod.DispersionOptimizer import DispersionOptimizer
from octopod.Cropper import Cropper
from octopod.Model import Model
from octopod.BScanAligner import BScanAligner
from octopod.Reporter import Reporter
from octopod.Flipper import Flipper

class Dataset:

    def __init__(self,raw_data_filename):
        """Initiate a Dataset object based on a .unp raw OCT data file. Pass just
        the base name if a Dataset is being created without an associated .unp file."""
        self.raw_data_filename = raw_data_filename
        self.logger = logging.getLogger(__name__)
        self.logger.info('Creating Dataset based on %s.'%raw_data_filename)
        raw_temp,raw_extension = os.path.splitext(self.raw_data_filename)
        self.hive_fn = raw_temp
        if self.hive_fn==self.raw_data_filename:
            sys.exit('Error: cannot create Hive with same name as raw data set.')
        self.xml_fn = raw_temp + '.xml'
        if not os.path.exists(self.xml_fn):
            xml_temp = raw_temp
            xml_temp = '_'.join(xml_temp.split('_')[:-1])+'.xml'
            if os.path.exists(xml_temp):
                self.logger.info('Could not locate %s.\nUsing %s instead.'%(self.xml_fn,xml_temp))
                self.xml_fn = xml_temp
        self.hive = Hive(self.hive_fn)


    def add_slo_frames(self,flist,width=512,height=512,dtype=np.uint16,skipbytes=0,labels=['SLO'],scaling_factor=1,invert=False):
        """Create a Dataset object from a series of SLO frames.
        This method is modeled slightly after Dataset.initialize, below."""
        
        self.hive.put('/config/n_vol',len(flist))
        self.hive.put('/config/n_slow',height)
        self.hive.put('/config/n_fast',width)
        self.hive.put('/config/n_depth',1)

        
        n_vol = self.hive.get('/config/n_vol')[0]
        n_slow = self.hive.get('/config/n_slow')[0]
        n_fast = self.hive.get('/config/n_fast')[0]

        raw_stores = {}
        for label in labels:
            raw_stores[label] = np.zeros((n_vol,n_slow,n_fast),dtype=dtype)

        for vol_index,f in enumerate(flist):
            print '%d of %d'%(vol_index+1,n_vol)
            with open(f,'r') as fid:
                fid.seek(skipbytes)
                arr = np.fromfile(fid,dtype=dtype)
                single_length = n_slow*n_fast
                full_length = single_length*len(labels)
                assert len(arr)==full_length
                for idx,label in enumerate(labels):
                    x1 = idx*single_length
                    x2 = x1+single_length
                    im = np.reshape(arr[x1:x2],(n_slow,n_fast))
                    raw_stores[label][vol_index,:,:] = im//scaling_factor
                    
            
        for label in labels:
            rs = raw_stores[label]
            if invert:
                rs = rs.max()-rs
            self.hive.put('projections/%s'%label,rs)
            

    def get_hive_handle(self):
        return self.hive

    def ecc_to_hive(self):
        eg = EccentricityGuesser()
        si_ecc,nt_ecc = eg.guess(self.raw_data_filename)
        self.hive.delete('eccentricity')
        self.hive.put('eccentricity/superior_inferior',si_ecc)
        self.hive.put('eccentricity/nasal_temporal',nt_ecc)
        self.hive.put('eccentricity/superior_and_nasal_are_negative',[np.nan])
        

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
            
        #idg.create_ids(self.hive,sn,at)

    def process(self):
        op = OCTProcessor(self.hive)
        op.run()

    def optimize_dispersion(self,n_lines=200,do_plot=False):
        do = DispersionOptimizer(self.hive)
        do.optimize(n_lines=n_lines,do_plot=do_plot)

    def flip(self):
        f = Flipper(self.hive)
        f.flip()

    def show(self):
        vols = self.hive.get('/processed_data')[:]
        nv,ns,nd,nf = vols.shape
        subvol = vols[0,:10,:,:]
        test = np.mean(np.abs(subvol),axis=0)
        plt.figure()
        plt.subplot(1,3,1)
        plt.imshow(test)
        plt.subplot(1,3,2)
        plt.imshow(np.log(test))
        plt.subplot(1,3,3)
        plt.plot(np.mean(test,axis=1))
        plt.show()
        

    def crop(self):
        c = Cropper(self.hive)
        c.crop()

    def model(self,keys=['ELM','ISOS','COST','RPE']):
        m = Model(self.hive)
        m.click_label(keys=keys)

    def align(self,do_plot=False):
        bsa = BScanAligner(self.hive)
        bsa.align_volumes(do_plot=do_plot)


    def label(self):
        m = Model(self.hive,data_block='flattened_data')
        m.write_axial_alignment()
        m.write_volume_labels()


    def report(self,do_log=False):
        report_directory = self.hive.filename.replace('.hdf5','')+'_report'
        r = Reporter(self.hive,report_directory=report_directory)
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

        # write parameters from the XML file to the hive file
        apf = AcquisitionParameterFile()
        apf.translate_xml_to_hive(self.xml_fn,self.hive)
        
        n_vol = self.hive.get('/config/n_vol')[()]
        n_slow = self.hive.get('/config/n_slow')[()]
        n_fast = self.hive.get('/config/n_fast')[()]
        n_depth = self.hive.get('/config/n_depth')[()]

        #raw_store = self.hive.make('raw_data',(n_vol,n_slow,n_fast,n_depth),dtype='u2')
        raw_store = np.zeros((n_vol,n_slow,n_fast,n_depth),dtype='u2')
        
        with open(self.raw_data_filename,'rb') as fid:
            for vol_index in range(n_vol):
                position = vol_index * n_depth * n_fast * n_slow * ocfg.raw_bytes_per_pixel
                fid.seek(position,0)
                vol = np.fromfile(fid,dtype=np.uint16,count=n_slow*n_depth*n_fast)
                vol = vol.reshape(n_slow,n_fast,n_depth)
                raw_store[vol_index,:,:,:] = vol
        self.hive['raw_data'] = raw_store
        
        p = ocfg.source_spectra[system_label] # polynomial for generating dataset's pixel->lambda function
        
        self.L = np.polyval(p,np.arange(n_depth))
        self.k_in = (2.0*np.pi)/self.L
        self.k_out = np.linspace(self.k_in[0],self.k_in[-1],n_depth)

        self.hive.put('L',self.L)
        self.hive.put('k_in',self.k_in)
        self.hive.put('k_out',self.k_out)
        self.hive.put('system_label',system_label)

        #self.hive.put('dispersion/coefficients',[0.0,0.0])
        
        self.ecc_to_hive()
        self.make_ids()
                
def test():
    ds = Dataset('./oct_test_volume/oct_test_volume_2T.unp')
    ds.initialize('2g_aooct')
        
if __name__=='__main__':
    test()
