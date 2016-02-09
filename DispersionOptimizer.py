import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from scipy import interpolate,optimize
import octopod_config as ocfg
import logging
from octopod.Processor import process

class DispersionOptimizer:

    def __init__(self,h5):
        self.h5 = h5
        self.logger = logging.getLogger(__name__)
        self.logger.info('Creating DispersionOptimizer object.')

    def make_test_frame(self,size=100):
        self.logger.info('Generating a random test frame from raw data.')
        n_depth = self.h5['/config/n_depth'][()]
        n_fast = self.h5['/config/n_fast'][()]
        n_slow = self.h5['/config/n_slow'][()]
        n_vol = self.h5['/config/n_vol'][()]

        test_frame = np.zeros((size,n_depth))
        hypercube = np.zeros(self.h5['raw_data'].shape)
        hypercube[...] = self.h5['raw_data']
        for k in range(size):
            v = np.random.randint(n_vol)
            s = np.random.randint(n_slow)
            f = np.random.randint(n_fast)
            test_frame[k,:] = hypercube[v][s][f][:]
        del hypercube
        return test_frame
        
        
    def process_frame(self,frame,c):
        test_frame = frame.copy()
        k_in = self.h5['k_in']
        k_out = self.h5['k_out']
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

    def optimize(self,test_frame):
        obj = lambda c_sub: self.dispersion_objective(test_frame,c_sub)
        c_sub0 = [0.0,0.0]
        bounds3 = [c_sub0[0]-1e-16,c_sub0[0]+1e-16]
        bounds2 = [c_sub0[1]-1e-10,c_sub0[1]+1e-10]
        
        result = optimize.brute(obj,(bounds3,bounds2),Ns=11,finish=None)
        
        bounds3a = (result[0]-1e-17,result[0]+1e-17)
        bounds2a = (result[1]-1e-11,result[1]+1e-11)
        result = optimize.brute(obj,(bounds3a,bounds2a),Ns=11,finish=None)
        
        c = [result[0],result[1],0.0,0.0]
        objective_value = obj(result)
        
        self.h5.require_group('dispersion')
        self.h5overwrite('dispersion/coefficients',np.array(c))
        self.h5overwrite('dispersion/objective_value',objective_value)
        self.h5overwrite('dispersion/bounds3',np.array(bounds3))
        self.h5overwrite('dispersion/bounds2',np.array(bounds2))
        
        return c
    
    def h5overwrite(self,key,value):
        try:
            del self.h5[key]
        except Exception as e:
            pass
        self.h5[key] = value
