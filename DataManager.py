import h5py
import hashlib

class IDGenerator:

    def create_ids(self,h5filename,subject_name,acquisition_date):
        # we need three ids in each hdf5 file: subject id, experiment id, and dataset id
        # we'll generate these by hashing the subject's name in Last_First format, a combination
        # of the name and date in YYYYMMDD format, and the first b-scan in the first volume of this dataset,
        # respectively

        def string_to_id(input_string,max_val=2**32):
            md5 = hashlib.md5()
            md5.update(input_string)
            return int(md5.hexdigest(),16)%max_val
        
        subject_id = string_to_id(subject_name)
        experiment_id = string_to_id(subject_name+acquisition_date)
        
        h5 = h5py.File(h5filename)

        test_frame = h5['raw_data'][0,0,:,:]
        dataset_id = string_to_id(test_frame)

        try:
            del h5['IDs']
        except Exception as e:
            pass

        h5.require_group('IDs')
        
