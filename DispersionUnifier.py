import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from scipy import interpolate,optimize
import octopod_config as ocfg
import logging,sys
from octopod.Processor import process
from octopod.DataStore import H5

class DispersionUnifier:

    def __init__(self,unp_list,auto=False):
        self.h5_list = []
        self.coefs = []
        self.k_ins = []
        self.k_outs = []
        self.n_slows = []

        winners = []
        
        for u in unp_list:
            h5fn = u.replace('.unp','')+'.hdf5'
            h5 = H5(h5fn)
            self.h5_list.append(h5)
            self.k_ins.append(h5['k_in'][:])
            self.k_outs.append(h5['k_out'][:])
            self.coefs.append(h5['dispersion/coefficients'][:])
            self.n_slows.append(h5['config/n_slow'].value)
            
        for h5,kin,kout,n_slow in zip(self.h5_list,self.k_ins,self.k_outs,self.n_slows):
            scores = []
            for idx,coef in enumerate(self.coefs):
                f = h5['raw_data'][0,int(n_slow/2),5:-5:,:]
                p = np.abs(process(f,kin,kout,coef))
                profile = np.mean(p,axis=1)
                pmax = np.max(profile)
                pmin = np.min(profile)
                pcontrast = (pmax - pmin)/(pmax + pmin)
                pgrad = np.max(np.abs(np.diff(profile)))
                c = 2
                score = pgrad*(pmax+c*pcontrast)
                scores.append(score)
                if not auto:
                    print '%d\t%0.1f\t%0.3f\t%0.1f\t%0.1f'%(idx+1,pmax,pcontrast,pgrad,score)
                    plt.figure()
                    plt.subplot(1,2,1)
                    clim = np.percentile(p,(5,99.9))
                    plt.imshow(p,interpolation='none',aspect='auto',clim=clim,cmap='gray')
                    plt.subplot(1,2,2)
                    plt.plot(profile)
            if not auto:
                plt.show()
                winner = int(raw_input('Choose best figure number: '))-1
            else:
                #print 'scores :',scores
                winner = np.argmax(scores)
                print 'winner: %d'%winner,self.coefs[winner]
            winners.append(winner)
        print 'winners:',winners
        mfw = max(set(winners), key=winners.count)
        print 'most frequent winner: ',mfw
