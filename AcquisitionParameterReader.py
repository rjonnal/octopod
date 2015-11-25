import octtools_config as ocfg
import h5py,os,sys
from bs4 import BeautifulSoup
from datetime import datetime
import logging

logging.basicConfig(level=logging.INFO)


class AcquisitionParameterReader:

    
    xml_dict = {}

    # populate xml_dict with required parameters from Yifan's XML grammar
    xml_dict['data_acquired_at'] = ('time_stamp',str)
    xml_dict['width'] = ('n_depth',int)
    xml_dict['height'] = ('n_fast',int)
    xml_dict['number_of_frames'] = ('n_slow',int)
    xml_dict['number_of_volumes'] = ('n_vol',int)
    xml_dict['x_scan_range'] = ('x_scan_mv',int)
    xml_dict['x_scan_offset'] = ('x_offset_mv',int)
    xml_dict['y_scan_range'] = ('y_scan_mv',int)
    xml_dict['y_scan_offset'] = ('y_offset_mv',int)
    xml_dict['number_of_bm_scans'] = ('n_bm_scans',int)
    xml_dict['c2'] = ('dispersion_c2',float)
    xml_dict['c3'] = ('dispersion_c3',float)
    
    def __init__(self):
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.info('Creating AcquisitionParameterReader.')

    
    def translate_xml_to_h5(self,fn,h5):
        self.logger.info('Creating "config" group in h5 file.')
        h5.create_group('config')
        h5['config'].create_dataset('filename',fn)

        fid = open(fn,'rb')
        soup = BeautifulSoup(fid.read(),"html.parser")
        headings = ['time','scanning_parameters','dispersion_parameters']
        for heading in headings:
            xml_list = soup.findAll(heading)
            for xml in xml_list:
                keys = xml.attrs.keys()
                for key in keys:
                    try:
                        h5key,op = self.xml_dict[key]
                        h5val = op(xml.attrs[key])
                        h5['config'].create_dataset(h5key,data=h5val)
                        self.logger.info('Setting variable %s to %s.'%(h5key,h5val))
                    except Exception as e:
                        self.logger.error(e)
                        sys.exit('Unrecoverable error: %s'%e)
                        
if __name__=='__main__':
    apfn = os.path.join('testing','2015-11-17-16-21-01-RE_3TR_0SR_20def_1.xml')
    f = h5py.File('./testing/foo.hdf5','w')
    apr = AcquisitionParameterReader()
    apr.translate_xml_to_h5(apfn,f)

