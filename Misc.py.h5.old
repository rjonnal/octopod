from octopod.DataStore import H5
import hashlib,os
import numpy as np

class EccentricityGuesser:
    '''Try to guess the eccentricity from the filename. Here are the assumptions:
    1. Abbreviations for temporal, nasal, superior, and inferior are upper-case first
    letters T, N, S, and I.
    2. If one of these letters is preceded by a string of characters that can be interpreted
    as a float, it's taken to mean the number of degrees in that direction.'''
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

    def guess_to_h5(self,h5filename):

        si_ecc,nt_ecc = self.guess(h5filename)
        
        h5 = H5(h5filename)
        h5.require_group('eccentricity')
        h5.put('eccentricity/superior_inferior',si_ecc)
        h5.put('eccentricity/nasal_temporal',nt_ecc)
        h5.put('eccentricity/superior_and_nasal_are_negative',[np.nan])
    
    def search_for_eccentricity(self,token,index):
        '''Try to convert characters just before index into a number.'''
        out = None
        for start in range(0,index):
            test = token[start:index]
            try:
                out = float(test)
                break
            except Exception as e:
                continue
        return out


    def token_guesser(self,token,letter):
        indices = self.findall(token,letter)
        value = None
        for idx in indices:
            value = self.search_for_eccentricity(token,idx)
            if value is not None:
                break
        return value
    
    def guess(self,filename):
        head,tail = os.path.split(filename)
        tokens = tail.split('_')
        s_ecc = []
        i_ecc = []
        n_ecc = []
        t_ecc = []
        for token in tokens:
            # find all instances of s, i, n, t
            s_ecc.append(self.token_guesser(token,'S'))
            i_ecc.append(self.token_guesser(token,'I'))
            n_ecc.append(self.token_guesser(token,'N'))
            t_ecc.append(self.token_guesser(token,'T'))

        # defaults in case no ecc is found
        si_ecc_out = 0.0
        nt_ecc_out = 0.0

        si_done = False
        for s in s_ecc:
            if s is not None:
                si_ecc_out = -s
                si_done = True
                break
                
        if not si_done:
            for i in i_ecc:
                if i is not None:
                    si_ecc_out = i
                    si_done = True
                    break

        nt_done = False
        for n in n_ecc:
            if n is not None:
                nt_ecc_out = -n
                nt_done = True
                break
                
        if not nt_done:
            for t in t_ecc:
                if t is not None:
                    nt_ecc_out = t
                    nt_done = True
                    break
        return si_ecc_out,nt_ecc_out

class IDGenerator:

    def create_ids(self,h5,subject_name,acquisition_date):
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
        
        test_frame = h5.get('raw_data')[0,0,:,:]
        dataset_id = string_to_id(test_frame)

        h5.require_group('IDs')
        h5.put('IDs/subject_id',subject_id)
        h5.put('IDs/experiment_id',experiment_id)
        h5.put('IDs/dataset_id',dataset_id)
        h5.put('acquisition_date',int(acquisition_date))

        #for key in h5['IDs'].keys():
            #print key,h5['IDs'][key].value
        #print h5['acquisition_date'].value


def test():
    h5filename = './oct_test_volume/oct_test_volume_2T.hdf5'
    eg = EccentricityGuesser()
    print 'Eccentricity guess: ',eg.guess(h5filename)
    eg.guess_to_h5(h5filename)
    h5 = H5(h5filename)

    idg = IDGenerator()
    idg.create_ids(h5,'jones_stick','20160202')
    
        
if __name__=='__main__':
    test()

