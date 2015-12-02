import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from scipy import interpolate,optimize
import octopod_config as ocfg
import logging

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
        test_frame = test_frame - np.mean(test_frame,axis=0)
        test_frame = test_frame.T
        k_in = self.h5['k_in']
        k_out = self.h5['k_out']
        k_interpolator = interpolate.interp1d(k_in,test_frame,axis=0,copy=False)
        test_frame = k_interpolator(k_out)
        
        dispersion_axis = k_out - np.mean(k_out)
        phase = np.exp(1j*np.polyval(c,dispersion_axis))
        test_frame = test_frame * phase[None].T
        test_frame = np.fft.fftshift(np.fft.fft(test_frame,axis=0),axes=0)

        n_depth = self.h5['/config/n_depth'][()]
        test_frame = test_frame[:n_depth/2,:]
        cutoff = ocfg.dc_cutoff
        test_frame = test_frame[:-cutoff,:]
        
        return test_frame
        
    def dispersion_objective(self,test_frame,c_sub):
        print c_sub,
        c_all = [c_sub[0],c_sub[1],0.0,0.0]
        frame = np.abs(self.process_frame(test_frame,c_all))
        out =  1.0/np.max(frame**2)
        print out
        return out

    def get(self,test_frame):
        obj = lambda c_sub: self.dispersion_objective(test_frame,c_sub)
        c_sub0 = [0.0,0.0]
        bounds3 = [c_sub0[0]-2e-16,c_sub0[0]+2e-16]
        bounds2 = [c_sub0[1]-2e-10,c_sub0[1]+2e-10]
        
        result = optimize.brute(obj,(bounds3,bounds2),Ns=41,finish=None)
        
        
        bounds3a = (result[0]-1e-17,result[0]+1e-17)
        bounds2a = (result[1]-1e-11,result[1]+1e-11)
        result = optimize.brute(obj,(bounds3a,bounds2a),Ns=11,finish=None)
        
        c = [result[0],result[1],0.0,0.0]
        return c
    #         elif method=='brute':
#             posttemp = result
#         else:
#             sys.exit('Please provide an optimization method.')

#         postcoef[0] = posttemp[0]
#         postcoef[1] = posttemp[1]

#         self.plotOptLog()
#         self.log.log('optimization time elapsed: %f'%t_elapsed)
#         self.optlog.log('# optimization time elapsed: %f'%t_elapsed)
#         self.showDispersionCompensation(coefs=postcoef,showPlot=False)
#         return postcoef
    
