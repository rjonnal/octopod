import h5py,os,sys
from bs4 import BeautifulSoup
from datetime import datetime
import logging

class AcquisitionParameterFile:
    
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
        self.logger = logging.getLogger(__name__)
        self.logger.info('Creating AcquisitionParameterFile object.')

    
    def translate_xml_to_h5(self,fn,h5):
        self.logger.info('Creating "config" group in h5 file.')
        h5.require_group('config')

        h5['config'].attrs['filename'] = fn

        fid = open(fn,'rb')
        soup = BeautifulSoup(fid.read(),"html.parser")
        headings = ['time','volume_size','scanning_parameters','dispersion_parameters']
        for heading in headings:
            xml_list = soup.findAll(heading)
            for xml in xml_list:
                keys = xml.attrs.keys()
                for key in keys:
                    try:
                        h5key,op = self.xml_dict[key]
                        h5val = op(xml.attrs[key])
                        try:
                            del h5['config'][h5key]
                        except Exception as e:
                            pass
                        h5['config'].create_dataset(h5key,data=h5val)
                        self.logger.info('Setting variable %s to %s.'%(h5key,h5val))
                    except Exception as e:
                        self.logger.error(e)
                        sys.exit('Unrecoverable error: %s'%e)


    def params_to_xml(self,category_dict):
        """category_dict is a dictionary of dictionaries. Keys should be the headings
        'time','volume_size','scanning_parameters','dispersion_parameters', and the values should
        be dictionaries containing key:value pairs to be written as attributes of the
        category's XML tag."""
        head = '<?xml version="1.0" encoding="utf-8"?>\n<MonsterList>\n\t<Monster>\n'
        tail = '\n\t<\Monster>\n<\MonsterList>'
        guts = ''
        for L1_key in category_dict.keys():
            attr_dict = category_dict[L1_key]
            guts = guts + '\n\t\t<%s '%L1_key
            for L2_key in attr_dict.keys():
                guts = guts + '%s="%s" '%(L2_key,attr_dict[L2_key])
            guts = guts + '/>'
        out = head+guts+tail
        self.logger.info('Generating XML: %s'%out.replace('\n',' '))
        return out

    def make_xml_file(self,data_filename,n_depth=2048,n_fast=1000,n_slow=100,n_vol=1):
        """This method generates Yifan-style XML files for bioptigen data sets."""
        timestamp = datetime.fromtimestamp(os.path.getmtime(data_filename)).strftime('%m/%d/%Y %I:%M:%S %p')
        head,tail = os.path.splitext(data_filename)
        out_filename = head + '.xml'
        if os.path.exists(out_filename):
            self.logger.error('Cannot overwrite %s. Please delete manually. Exiting.'%out_filename)
            sys.exit()
        test = {}
        test['time']={'Data_Acquired_at':timestamp}
        test['volume_size']={'width':'%d'%n_depth,'height':'%d'%n_fast,'number_of_frames':'%d'%n_slow,'number_of_volumes':'%d'%n_vol}
        test['scanning_parameters']={'x_scan_range':'3168','x_scan_offset':'0','y_scan_range':'2112','y_scan_offset':'0','number_of_bm_scans':'1'}
        test['dispersion_parameters']={'c2':'0.0','c3':'0.0'}
        # perform a sanity check on file size
        bytecount = os.stat(data_filename).st_size
        expected_bytecount = int(test['volume_size']['width'])*int(test['volume_size']['height'])*int(test['volume_size']['number_of_frames'])*int(test['volume_size']['number_of_volumes'])*2
        if not bytecount==expected_bytecount:
            self.logger.error('Computed bytecount and expected bytecount for %s do not agree. Exiting.'%data_filename)
            sys.exit()
        with open(out_filename,'w') as fid:
            fid.write(self.params_to_xml(test))
        
            
            
    
if __name__=='__main__':
    apfn = os.path.join('testing','2015-11-17-16-21-01-RE_3TR_0SR_20def_1.xml')
    f = h5py.File('./testing/foo.hdf5','w')
    apr = AcquisitionParameterReader()
    apr.translate_xml_to_h5(apfn,f)

