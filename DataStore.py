'''A thin wrapper for h5py, in order to abstract away from the latter's implementation details. This class is meant to
guard the bulk of the code from the data archiving implementation, such that if, for instance, we decided to switch from
HDF5 to SQL, we could implement it quickly.
'''

import h5py
import logging
logging.basicConfig(level='INFO')

class DataStore:

    '''An abstract interface for data storage. The assumptions: 

    1. the storage system is heirarchical, with branches called groups (analgous to folders or directories) and
    leaves called datasets (analogous to files).

    2. the data to put into storage are either numpy arrays (or equivalent structures such as Python lists of floats or
    ints) or strings.

    '''

    def __init__(self,filename):
        self.filename = filename
        self.logger = logging.getLogger(__name__)
        self.logger.info('Creating %s object.'%__name__)

    def put(self,location,data,short_descriptor=None):
        '''Put some data into the store.

        Args: 

            location (str): the location for the data (may be heirarchical or not,
            e.g. '/not/very/important/data/juggling_records' or '/blood_pressure')

            data (Numpy array or str): the data to store

            short_descriptor (str): a short description of the data
        '''
        pass

    def get(self,location):
        '''Get some data out of the store.

        Args:

            location (str): the location of the sought data

        Returns:

            output (various types): the referenced data, either Numpy array or str
        '''
        pass

    def __str__(self):
        return 'DataStore object located at %s.'%self.filename


    def make(self,dims):
        '''Make an empty data structure with dimensions specified by the tuple
        dims.

        Args:

            dims (tuple of ints): the dimensions of the desired empty Numpy array

        Returns:

            output (Numpy array): an interface into the desired array
        '''


    def close(self):
        pass

class H5(DataStore):
    
    def __init__(self,filename):
        DataStore.__init__(self,filename)
        self.h5 = h5py(filename)

    def put(self,location,data,short_descriptor=None):
        try:
            del self.h5[location]
        except Exception as e:
            pass

        self.h5.create_dataset(location,data=data)


    def make(self,dims,dtype=''):
        

    
