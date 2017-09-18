'''A thin wrapper for py, in order to abstract away from the latter's implementation details. This class is meant to
guard the bulk of the code from the data archiving implementation, such that if, for instance, we decided to switch from
HDF5 to SQL, we could implement it quickly.
'''

import h5py
import os,sys,glob
import shutil
import logging
import numpy as np
logging.basicConfig(level='INFO')

class DataStore:

    '''An abstract interface for data storage. The assumptions: 

    1. the storage system is heirarchical, with branches called groups (analgous to folders or directories) and
    leaves called datasets (analogous to files).

    2. the data to put into storage are either numpy arrays (or equivalent structures such as Python lists of floats or
    ints) or strings.

    '''

    def __init__(self,filename,mode='w'):
        self.filename = filename
        self.logger = logging.getLogger(__name__)
        self.logger.info('Creating %s object from %s.'%(self.__class__,filename))
        self.mode = mode

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
    
    def __init__(self,filename,mode='w'):
        DataStore.__init__(self,filename,mode)
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
        if self.mode.find('w')==-1:
            self.logger.info('Mode is %s; cannot repack. Exiting.'%self.mode)
            sys.exit()
        newh5fn = self.filename+'.repacked.hdf5'
        if os.path.exists(newh5fn):
            self.move(newh5fn,newh5fn+'.garbage')
            
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

    def print_helper(self,thing,depth=0,max_depth=np.inf):
        if depth>=max_depth:
            return
        try:
            print thing.keys()
            for key in thing.keys():
                try:
                    print '\t'*depth,key,':',
                    self.print_helper(thing[key],depth+1,max_depth)
                    #print
                except:
                    pass
        except:
            print '(leaf)',thing.shape

    def catalog(self,max_depth=np.inf):
        self.print_helper(self.h5,max_depth=max_depth)

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
        if self.mode.find('w')==-1:
            self.logger.info('Mode is %s; cannot put. Exiting.'%self.mode)
            sys.exit()
            
        try:
            del self.h5[location]
        except Exception as e:
            pass

        self.h5.create_dataset(location,data=data)
        self.logger.info('Putting %s into H5 file at %s.'%(self.h5[location],location))
        self.write_descriptor(location,short_descriptor)


    def require_group(self,location):
        if self.mode.find('w')==-1:
            self.logger.info('Mode is %s; cannot require group. Exiting.'%self.mode)
            sys.exit()
        self.h5.require_group(location)


    def get_shape(self,location):
        return self.h5[location].shape
    
    def make(self,location,dims,dtype='f8',short_descriptor=None):
        #"<i1", "<i2", "<i4", "<i8", ">i1", ">i2", ">i4", ">i8", "|i1", "|u1", 
        #"<u1", "<u2", "<u4", "<u8", ">u1", ">u2", ">u4", ">u8",
        #"<f4", "<f8", ">f4", ">f8", "<c8", "<c16", ">c8", ">c16",
        #"|S1", "|S2", "|S33", "|V1", "|V2", "|V33"]
        if self.mode.find('w')==-1:
            self.logger.info('Mode is %s; cannot make. Exiting.'%self.mode)
            sys.exit()
            
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
        if self.mode.find('w')==-1:
            self.logger.info('Mode is %s; cannot write descriptor. Exiting.'%self.mode)
            sys.exit()
        if descriptor is not None:
            location = '%s_IS_%s'%(location,descriptor.replace(' ','_'))
            self.put(location,[1])

    def write_attribute(self,group,attr_key,attr_value):
        if self.mode.find('w')==-1:
            self.logger.info('Mode is %s; cannot write attribute. Exiting.'%self.mode)
            sys.exit()
        self.h5.require_group(group)
        self.h5[group].attrs[attr_key] = attr_value

    def delete(self,key):
        if self.mode.find('w')==-1:
            self.logger.info('Mode is %s; cannot delete. Exiting.'%self.mode)
            sys.exit()
        try:
            del self.h5[key]
        except Exception as e:
            pass

class Hive(DataStore):
    
    def __init__(self,root_location,mode='w'):
        self.logger = logging.getLogger(__name__)
        self.root_location = root_location
        self.mode = mode
        try:
            os.makedirs(self.root_location)
            self.logger.info('Creating Hive at %s.'%self.root_location)
        except Exception as e:
            if os.path.exists(self.root_location):
                self.logger.info('Using existing Hive at %s.'%self.root_location)
            else:
                sys.exit(e)
        
    def move(self,src,dest):
        self.logger.info('move: moving %s -> %s'%(src,dest))
        try:
            shutil.move(src,dest)
        except Exception as e:
            self.logger.error(e)
            sys.exit(e)

    def repack(self):
        pass

    def keys(self):
        return glob.glob(os.path.join(self.root_location,'*'))

    def print_helper(self,thing,depth=0,max_depth=np.inf):
        if depth>=max_depth:
            return
        try:
            thingkeys = glob.glob(os.path.join(thing,'*'))
            print thingkeys
            for key in thingkeys:
                try:
                    print '\t'*depth,key,':',
                    self.print_helper(os.path.join(os.path.join(thing,key),'*'),depth+1,max_depth)
                except:
                    pass
        except:
            print '(leaf)',thing

    def catalog(self,max_depth=np.inf):
        print 'Catalog:'
        self.print_helper(self.root_location,max_depth=max_depth)
        print
        
    def has(self,key):
        return os.path.exists(os.path.join(self.root_location,key+'.npy'))
    
    def __getitem__(self,key):
        if key[0]=='/':
            key = key[1:]
        fn = os.path.join(self.root_location,key)
        if os.path.exists(fn+'.npy'):
            return np.load(fn+'.npy')
        else:
            return Hive(fn)

    def put(self,location,data,short_descriptor=None,overwrite=False):
        if location[0]=='/':
            location = location[1:]

        if type(data)==list:
            data = np.array(data)
        elif not type(data)==np.ndarray:
            data = np.array([data])
            
        if self.mode.find('w')==-1:
            self.logger.info('Mode is %s; cannot put. Exiting.'%self.mode)
            sys.exit()
            
        fn = os.path.join(self.root_location,location)+'.npy'

        if os.path.exists(fn) and not overwrite:
            self.logger.info('%s already exists; to overwrite call Hive.put with overwrite=True.'%location)
            return

        path,shortfn = os.path.split(fn)
        if not os.path.exists(path):
            os.makedirs(path)
        
        try:
            os.remove(fn)
        except Exception as e:
            self.logger.info(e)

        np.save(fn,data)
        if short_descriptor is not None:
            txtfn = os.path.join(self.root_location,location)+'.txt'
            with open(txtfn,'wb') as fid:
                fid.write(short_descriptor)
        self.logger.info('Putting %s (%s,%s) into H5 file at %s.'%(type(data),list(data.shape),data.dtype,location))


    def get_shape(self,location):
        return self[location].shape
    
    def get(self,location):
        return self[location]
    
    def close(self):
        pass
    
    def delete(self,key):
        if self.mode.find('w')==-1:
            self.logger.info('Mode is %s; cannot delete. Exiting.'%self.mode)
            sys.exit()
        fn = os.path.join(self.root_location,key)+'.npy'
        try:
            os.remove(fn)
        except Exception as e:
            pass
