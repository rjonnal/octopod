import numpy as np
import scipy as sp
from matplotlib import pyplot as plt

class DispersionOptimizer:

    def __init__(self,initial_coefficients=[0.0,0.0,0.0,0.0]):
        self.initial_coefficients = initial_coefficients


    def dispersion_objective(self,test_frame,c):
        
    def get(self,test_frame):
        test_frame = self.makeTestFrame(N,width)
        
        obj = lambda c: self.dispersionObjective(test_frame,c,plot=plot)

#        self.log.log('Getting initial values for optimization from...')
#        c0 = self.readGlobalDispersionCoefs(self.isNfl)
#        if c0 is None:
#            c0 = self.dispersionCoefs[:2]

        c0 = [0.0,0.0]

        precoef = [0.,0.,0.,0.]
        postcoef = [0.,0.,0.,0.]

        #lowers = ocfg.DISPERSION_COEFS_LBOUNDS[:2]
        #uppers = ocfg.DISPERSION_COEFS_UBOUNDS[:2]
        lowers = [c0[0]-2e-16,c0[1]-2e-10]
        uppers = [c0[0]+2e-16,c0[1]+2e-10]
        
        bounds3 = (lowers[0],uppers[0])
        bounds2 = (lowers[1],uppers[1])


        if plot:
            plt.figure()

        t0 = time()
        self.log.log('optimization start')
        self.optlog.log('# optimization start')
        self.optlog.log('# starting with coefs %0.3e,%0.3e'%(c0[0],c0[1]))
        self.optlog.log('# order3\torder2\t\tmetric')


        method='brute'

        if method=='brute':
            result = optimize.brute(obj,(bounds3,bounds2),Ns=41,finish=None)
            bounds3a = (result[0]-1e-17,result[0]+1e-17)
            bounds2a = (result[1]-1e-11,result[1]+1e-11)
            result = optimize.brute(obj,(bounds3a,bounds2a),Ns=11,finish=None)
        elif method=='anneal':
            result = optimize.anneal(obj,x0=c0,schedule='fast',lower=lowers,upper=uppers,maxeval=None)
        else:
            sys.exit('Please provide an optimization method.')
            
        t_elapsed = time() - t0

        if plot:
            plt.ioff()
            plt.close()

        if method=='anneal':
            posttemp = result[0]
        elif method=='brute':
            posttemp = result
        else:
            sys.exit('Please provide an optimization method.')

        postcoef[0] = posttemp[0]
        postcoef[1] = posttemp[1]

        self.plotOptLog()
        self.log.log('optimization time elapsed: %f'%t_elapsed)
        self.optlog.log('# optimization time elapsed: %f'%t_elapsed)
        self.showDispersionCompensation(coefs=postcoef,showPlot=False)
        return postcoef
    
