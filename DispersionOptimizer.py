import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from scipy import interpolate,optimize
import octopod_config as ocfg
import logging
from octopod.Processor import process
from octopod.DataStore import H5

class DispersionOptimizer:

    def __init__(self,h5):
        self.h5 = h5
        self.logger = logging.getLogger(__name__)
        self.logger.info('Creating DispersionOptimizer object.')

    def make_test_frame(self,size=100):
        self.logger.info('Generating a random test frame from raw data.')
        
        n_depth = self.h5.get('/config/n_depth')[()]
        n_fast = self.h5.get('/config/n_fast')[()]
        n_slow = self.h5.get('/config/n_slow')[()]
        n_vol = self.h5.get('/config/n_vol')[()]

        test_frame = np.zeros((size,n_depth))
        hypercube = np.zeros(self.h5.get('raw_data').shape)
        hypercube[...] = self.h5.get('raw_data')[...]
        for k in range(size):
            v = np.random.randint(n_vol)
            s = np.random.randint(n_slow)
            f = np.random.randint(n_fast)
            test_frame[k,:] = hypercube[v][s][f][:]
        del hypercube
        return test_frame
        
    def process_frame(self,frame,c):
        test_frame = frame.copy()
        k_in = self.h5.get('k_in')
        k_out = self.h5.get('k_out')
        test_frame = process(frame,k_in,k_out,c)
        return test_frame
        
    def dispersion_objective(self,test_frame,c_sub):
        print c_sub,
        c_all = [c_sub[0],c_sub[1],0.0,0.0]
        frame = np.abs(self.process_frame(test_frame,c_all))
        colmax = np.max(frame**2,axis=0)
        out =  1.0/np.mean(colmax)
        print out
        return out

    def optimize(self,test_frame=None,dry_run=False):
        if test_frame is None:
            test_frame = self.make_test_frame(1000)
            
        obj = lambda c_sub: self.dispersion_objective(test_frame,c_sub)
        c_sub0 = [0.0,0.0]
        bounds3 = [c_sub0[0]-1e-17,c_sub0[0]+1e-17]
        bounds2 = [c_sub0[1]-1e-11,c_sub0[1]+1e-11]
        
        result = optimize.brute(obj,(bounds3,bounds2),Ns=11,finish=None)
        
        bounds3a = (result[0]-1e-18,result[0]+1e-18)
        bounds2a = (result[1]-1e-12,result[1]+1e-12)
        result = optimize.brute(obj,(bounds3a,bounds2a),Ns=11,finish=None)
        
        c = [result[0],result[1],0.0,0.0]
        objective_value = obj(result)

        if not dry_run:
            self.h5.require_group('dispersion')
            self.h5.put('dispersion/coefficients',np.array(c))
            self.h5.put('dispersion/objective_value',objective_value)
            self.h5.put('dispersion/bounds3',np.array(bounds3))
            self.h5.put('dispersion/bounds2',np.array(bounds2))
        
        return c
    

if __name__=='__main__':

    do = DispersionOptimizer('./oct_test_volume/oct_test_volume_2T.hdf5')
    test_frame = do.make_test_frame()
    print do.optimize(test_frame)
