import os,sys
import fnmatch
import logging
import octopod_config as ocfg
from octopod import H5

logging.basicConfig(level=logging.DEBUG)

class FileManager:

    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.logger.info('Creating FileManager.')
        paths = ocfg.data_paths
        self.keys = paths.keys()

        self.raw_files = {}
        self.h5_files = {}

        
        for key in self.keys:
            self.raw_files[key] = []
            self.h5_files[key] = []
            for root,dirnames,filenames in os.walk(paths[key]) :
                for filename in fnmatch.filter(filenames,'*.%s'%(ocfg.raw_data_extension)):
                    full_fn = os.path.join(root,filename)
                    self.raw_files[key].append(full_fn)
                    #self.logger.info('Adding %s to set %s.'%(full_fn,key))
                    
                for filename in fnmatch.filter(filenames,'*.%s'%(ocfg.h5_extension)):
                    full_fn = os.path.join(root,filename)
                    self.h5_files[key].append(full_fn)
                    self.logger.info('Adding %s to set %s.'%(full_fn,key))


    def report(self,filename='./octopod_report.csv'):

        logging.basicConfig(level=logging.WARNING)
        
        headings = ['subject_name','date','time_stamp','subject_id','experiment_id','dataset_id',
                    'processed','fast_scan_deg','n_fast','fast_sampling_um','slow_scan_deg',
                    'n_slow','slow_sampling_um','dispersion_optimized','c3','c2','model_done','model_labels']

        formats = ['%s','%d','%s','%d','%d','%d',
                   '%s','%0.1f','%d','%0.2f','%0.1f',
                   '%d','%0.2f','%s','%0.3e','%0.3e','%s','%s']

        fid = open(filename,'w')
        fid.write(', '.join(headings)+',\n')
            
        
        for key in self.keys:
            for filename in self.h5_files[key]:
                row = {}
                h5 = H5(filename)
                hkeys = h5.keys()

                # try to get subject_name
                path,fn = os.path.split(filename)
                try:
                    snfid = open(os.path.join(path,'subject_name.txt'))
                    sn = snfid.read().strip()
                    snfid.close()
                except Exception as e:
                    print e
                    sn = 'unknown'

                row['subject_name'] = sn

                # try to get acquisition date
                try:
                    a_date = h5.get('acquisition_date').value
                except:
                    a_date = 'unknown'

                row['date'] = a_date

                # try to get IDs
                IDs = h5.get('IDs')
                try:
                    sid = IDs['subject_id'].value
                except:
                    sid = 'unknown'
                try:
                    did = IDs['dataset_id'].value
                except:
                    did = 'unknown'
                try:
                    eid = IDs['experiment_id'].value
                except:
                    eid = 'unknown'

                row['subject_id'] = sid
                row['dataset_id'] = did
                row['experiment_id'] = eid

                # make sure processed and size matches config numbers
                try:
                    nv,ns,nd,nf = h5.get('processed_data').shape
                    processed = True
                except:
                    nv,ns,nd,nf = -1,-1,-1,-1
                    processed = False

                try:
                    processed = processed and (h5.get('config')['n_slow'].value==ns and
                        h5.get('config')['n_vol'].value==nv and
                        h5.get('config')['n_fast'].value==nf and
                        h5.get('config')['n_depth'].value==nd*2)
                except:
                    processed = False

                row['processed'] = processed

                # get scanning parameters
                try:
                    x_scan_mv = h5.get('config')['x_scan_mv'].value
                    y_scan_mv = h5.get('config')['y_scan_mv'].value
                    x_scan = x_scan_mv/ocfg.x_mv_per_deg
                    y_scan = y_scan_mv/ocfg.y_mv_per_deg
                except:
                    x_scan = -1
                    y_scan = -1

                row['fast_scan_deg'] = x_scan
                row['slow_scan_deg'] = y_scan

                row['n_fast'] = nf
                row['n_slow'] = ns

                fsu = 300.0*x_scan/nf
                ssu = 300.0*y_scan/ns

                row['fast_sampling_um'] = fsu
                row['slow_sampling_um'] = ssu
                
                # get time stamp
                try:
                    time_stamp = h5.get('config')['time_stamp'].value
                except:
                    time_stamp = 'unknown'

                row['time_stamp'] = time_stamp
                
                # get dispersion optimization info
                try:
                    coefs = h5.get('dispersion')['coefficients'][:]
                    dispersion_done = True
                except:
                    dispersion_done = False
                    coefs = [-1,-1,-1,-1]

                row['dispersion_optimized'] = dispersion_done
                row['c3'] = coefs[0]
                row['c2'] = coefs[1]


                # get model info
                try:
                    model = h5.get('model')['profile'][:]
                    model_done = True
                except:
                    model_done = False

                row['model_done'] = model_done

                model_labels = ''
                if model_done:
                    try:
                        model_labels = ' '.join(h5.get('model')['labels'].keys())
                    except:
                        pass

                row['model_labels'] = model_labels


                for h,f in zip(headings,formats):
                    fid.write(('%s, '%f)%(row[h]))
                fid.write(',\n')

        fid.close()
                

                #print h5.get('model').keys()
                
                #print h5.get('dispersion').keys()
                #print hkeys
                #print h5.get('IDs').keys()

        
                    
    def get(self,key,terms=[]):
        try:
            files = self.files[key]
            if len(terms):
                out = []
                for f in files:
                    for term in terms:
                        if f.lower().find(term.lower())>-1:
                            out.append(f)
            else:
                out = files
            return out
                        
        except KeyError as e:
            self.logger.info('KeyError: "%s". Returning [].'%(e.message))
            return []

    def get_one(self,key,terms=[]):
        files = self.get(key)
        if len(files):
            matches = 0
            winner = 0
            for idx,f in enumerate(files):
                temp_matches = 0
                for term in terms:
                    if f.lower().find(term.lower())>-1:
                        temp_matches = temp_matches + 1
                if temp_matches>matches:
                    matches = temp_matches
                    winner = idx
            return files[winner]



if __name__=='__main__':

    fm = FileManager()
    fm.report()
