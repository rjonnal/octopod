import h5py
import hashlib,os

class EccentricityGuesser:

    def findall(self,target,pattern):
        out = []
        remainder = target
        location = 0
        while len(remainder)>0 and location>-1:
            location = remainder.find(pattern)
            remainder = remainder[location+1:]
            if location>-1:
                out.append(location)
        return out

    def search_for_eccentricity(self,token,index):
        '''Try to convert characters just before index into a number.'''
        out = -1
        for start in range(0,index-1):
            test = token[start:index]
            try:
                out = int(test)
                break
            except Exception as e:
                continue

        return out


    def token_guesser(self,token,letter):
        indices = self.findall(token.lower(),letter.lower())
        out = []
        for idx in indices:
            value = self.search_for_eccentricity(token,idx)
            out.append(value)
        print token,letter,out
        return out
    
    def guess(self,filename):
        head,tail = os.path.split(filename)
        tokens = tail.split('_')
        for token in tokens:
            # find all instances of s, i, n, t
            s_idx = self.token_guesser(token,'s')
            i_idx = self.token_guesser(token,'i')
            n_idx = self.token_guesser(token,'n')
            t_idx = self.token_guesser(token,'t')
            
            
        

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

        h5.create_group('IDs')
        h5['IDs'].create_dataset('subject_id',data=subject_id)
        h5['IDs'].create_dataset('experiment_id',data=experiment_id)
        h5['IDs'].create_dataset('dataset_id',data=dataset_id)

        try:
            del h5['acquisition_date']
        except Exception as e:
            pass

        h5.create_dataset('acquisition_date',data=int(acquisition_date))

        #for key in h5['IDs'].keys():
            #print key,h5['IDs'][key].value
        #print h5['acquisition_date'].value
