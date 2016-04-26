'''A thin wrapper for h5py, in order to abstract away from the latter's implementation details. This class is meant to
guard the bulk of the code from the data archiving implementation, such that if, for instance, we decided to switch from
HDF5 to SQL, we could implement it quickly.
'''

import h5py
import os
import shutil
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
        self.logger.info('Creating %s object from %s.'%(self.__class__,filename))


    def __getitem__(self,key):
        return self.get(key)

    
    def close(self):
        '''Close interface to this data store.'''
        pass
    
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


class H5(DataStore):
    
    def __init__(self,filename):
        DataStore.__init__(self,filename)
        self.h5 = h5py.File(filename)
        self.filename = filename

    def move(self,src,dest):
        self.logger.info('move: moving %s -> %s'%(src,dest))
        try:
            shutil.move(src,dest)
        except Exception as e:
            self.logger.error(e)
            sys.exit(e)

        
    def repack(self):
        newh5fn = self.filename+'.repacked.hdf5'
        if os.path.exists(newh5fn):
            sys.exit('repack: Please delete %s before calling H5.repack'%newh5fn)
            
        newh5 = h5py.File(newh5fn)
        for key in self.h5.keys():
            h5py.h5o.copy(self.h5.id,key,newh5.id,key)
            self.logger.info('repack: Copying %s from %s to %s.'%(key,self.filename,newh5fn))
        
        self.h5.flush()
        newh5.flush()
        newh5.close()
        self.h5.close()
        
        self.move(self.filename,self.filename+'.repack.backup')
        self.move(newh5fn,self.filename)
        
        self.h5 = h5py.File(self.filename)
        

    def keys(self):
        return self.h5.keys()

    def has(self,key):
        out = True
        try:
            junk = self.h5[key]
        except:
            out = False
        return out
    # While convenient, implementation of getitem greatly reduces the flexibility of
    # the DataStore object, since alternative implementations (e.g. SQl) won't have
    # dictionary behavior by default.
    #def __getitem__(self,key):
    #    return self.get(key)

    def put(self,location,data,short_descriptor=None):
        try:
            del self.h5[location]
        except Exception as e:
            pass

        self.h5.create_dataset(location,data=data)
        self.logger.info('Putting %s into H5 file at %s.'%(self.h5[location],location))
        self.write_descriptor(location,short_descriptor)


    def require_group(self,location):
        self.h5.require_group(location)


    def get_shape(self,location):
        return self.h5[location].shape
    
    def make(self,location,dims,dtype='f8',short_descriptor=None):
        #"<i1", "<i2", "<i4", "<i8", ">i1", ">i2", ">i4", ">i8", "|i1", "|u1", 
        #"<u1", "<u2", "<u4", "<u8", ">u1", ">u2", ">u4", ">u8",
        #"<f4", "<f8", ">f4", ">f8", "<c8", "<c16", ">c8", ">c16",
        #"|S1", "|S2", "|S33", "|V1", "|V2", "|V33"]
        try:
            del self.h5[location]
        except Exception as e:
            pass

        self.h5.create_dataset(location,dims,dtype=dtype)
        
        return self.h5[location]

    def get(self,location):
        return self.h5[location]
    
    def close(self):
        self.h5.close()
    
    def write_descriptor(self,location,descriptor):
        if descriptor is not None:
            location = '%s_IS_%s'%(location,descriptor.replace(' ','_'))
            self.put(location,[1])

    def write_attribute(self,group,attr_key,attr_value):
        self.h5.require_group(group)
        self.h5[group].attrs[attr_key] = attr_value

    def delete(self,key):
        try:
            del self.h5[key]
        except Exception as e:
            pass
