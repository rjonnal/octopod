import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from scipy import interpolate,optimize
import octopod_config as ocfg
import logging,sys
from octopod.Processor import process
from octopod.DataStore import H5

class DispersionOptimizer:

    def __init__(self,h5):
        if type(h5)==str:
            self.h5 = H5(h5)
        else:
            self.h5 = h5
        self.logger = logging.getLogger(__name__)
        self.logger.info('Creating DispersionOptimizer object.')
        self.Ns = 21

    def show_dispersion_results(self,N=None,save=False):
        if N is None:
            N = 2*self.Ns

        logs = []
        log_names = []
        for key in self.h5.get('dispersion/logs').keys():
            logs.append(self.h5.get('dispersion/logs/%s'%key))
            log_names.append(key)
        try:
            coefs = self.h5.get('dispersion/coefficients')
        except Exception as e:
            self.logger.info(e)

        def plot(LL,name=''):
            y = LL[:,0]
            x = LL[:,1]
            z = LL[:,2]
            ymax = np.max(y)
            ymin = np.min(y)
            xmax = np.max(x)
            xmin = np.min(x)
            XX,YY = np.mgrid[xmin:xmax:N*1j,ymin:ymax:N*1j]
            grid = interpolate.griddata(np.array([x,y]).T,z,(XX,YY),method='cubic')
            plt.pcolormesh(XX,YY,grid)
            plt.plot(coefs[1],coefs[0],'ws')
            plt.axis('tight')
            plt.colorbar()
            plt.title(name)

        nlog = len(logs)
        plt.figure(figsize=(4*nlog,4))
        for idx,(log,log_name) in enumerate(zip(logs,log_names)):
            plt.subplot(1,nlog,1+idx)
            plot(log,log_name)

        if save:
            outfn = self.h5.filename.replace('.hdf5','')+'_dispersion_optimization.png'
            plt.savefig(outfn)
            plt.pause(1)
            plt.close()
        else:
            plt.show()

    def make_test_frame(self,size=1000):
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
        test_frame = process(frame,k_in,k_out,c)[:-ocfg.dc_cutoff,:]
        return test_frame
        
    def dispersion_objective(self,test_frame,c_sub,log=[],do_plot=False):
        c_all = [c_sub[0],c_sub[1],0.0,0.0]
        frame = np.abs(self.process_frame(test_frame,c_all))
        colmax = np.max(frame**2,axis=0)
        out =  1.0/np.mean(colmax)
        
        last_pct_done = float(len(log))/float(self.Ns**2)
        log.append([c_sub[0],c_sub[1],out])
        pct_done = float(len(log))/float(self.Ns**2)
        if np.floor(pct_done*10)-np.floor(last_pct_done*10):
            self.logger.info('%d percent done.'%(pct_done*100))
        
        if do_plot:
            plt.subplot(1,2,1)
            plt.cla()
            plt.imshow(frame,interpolation='none',aspect='auto')
            plt.title('%0.1e'%out)
            plt.subplot(1,2,2)
            plt.cla()
            plt.hist(np.ravel(frame))
            plt.title('%0.1e,%0.1e'%(c_sub[0],c_sub[1]))
            plt.pause(.001)
            
        return out

    def optimize(self,test_frame=None,dry_run=False,do_plot=False,method='brute',n_lines=1000):
        if test_frame is None:
            test_frame = self.make_test_frame(n_lines)
            
        coarse_log = []
        fine_log = []
        
        obj = lambda c_sub: self.dispersion_objective(test_frame,c_sub,coarse_log,do_plot)
        c_sub0 = [0.0,0.0]
        bounds3 = [c_sub0[0]-1e-17,c_sub0[0]+1e-17]
        bounds2 = [c_sub0[1]-1e-11,c_sub0[1]+1e-11]

        
        step3 = (bounds3[1]-bounds3[0])/float(self.Ns)
        step2 = (bounds2[1]-bounds2[0])/float(self.Ns)
        self.logger.info('Starting coarse optimization.')
        self.logger.info('Bounds (3rd order): %0.1e,%0.1e'%tuple(bounds3))
        self.logger.info('Bounds (2nd order): %0.1e,%0.1e'%tuple(bounds2))
        result = optimize.brute(obj,(bounds3,bounds2),Ns=self.Ns,finish=None)
        self.logger.info('Optimum: %0.1e,%0.1e'%tuple(result))
        
        
        obj = lambda c_sub: self.dispersion_objective(test_frame,c_sub,fine_log,do_plot)
        bounds3a = (result[0]-step3,result[0]+step3)
        bounds2a = (result[1]-step2,result[1]+step2)
        self.logger.info('Starting fine optimization.')
        self.logger.info('Bounds (3rd order): %0.1e,%0.1e'%tuple(bounds3a))
        self.logger.info('Bounds (2nd order): %0.1e,%0.1e'%tuple(bounds2a))
        result = optimize.brute(obj,(bounds3a,bounds2a),Ns=self.Ns,finish=None)
        self.logger.info('Optimum: %0.1e,%0.1e'%tuple(result))
        
        c = [result[0],result[1],0.0,0.0]
        objective_value = obj(result)

        if not dry_run:
            self.h5.require_group('dispersion')
            self.h5.put('dispersion/coefficients',np.array(c))
            self.h5.put('dispersion/objective_value',objective_value)
            self.h5.put('dispersion/bounds3',np.array(bounds3))
            self.h5.put('dispersion/bounds2',np.array(bounds2))
            self.h5.put('dispersion/logs/coarse',np.array(coarse_log))
            self.h5.put('dispersion/logs/fine',np.array(fine_log))
        self.logger.info('optimized coefficients: %0.2e, %0.2e'%(result[0],result[1]))
        return c
    

if __name__=='__main__':

    do = DispersionOptimizer('./oct_test_volume/oct_test_volume_2T.hdf5')
    test_frame = do.make_test_frame(100)
    do.optimize(test_frame,do_plot=True)
    do.show_dispersion_results()
