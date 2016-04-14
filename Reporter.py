from octopod.DataStore import H5
import os,sys
from utils import get_z_sampling
import logging
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate
logging.basicConfig(level=logging.DEBUG)

class Reporter:

    def __init__(self,h5,report_directory='.'):
        self.logger = logging.getLogger(__name__)
        self.h5 = h5
        self.report_directory = report_directory
        self.makedirs(report_directory)


    def makedirs(self,directory_name):
        if not os.path.exists(directory_name):
            os.makedirs(directory_name)
    

    def model_report(self):
        tag = 'axial scattering model'
        try:
            working_profile = self.h5.get('model/profile')[:]
        except Exception as e:
            return

        try:
            L = self.h5.get('L')[:]
            dz = abs(get_z_sampling(L[0],L[-1]))*1e6
            #x = xbarfoo
        except Exception as e:
            dz = 1.0

        z = np.arange(len(working_profile))*dz

        label_dict = {}
        try:
            keys = self.h5.get('model/labels').keys()
        except Exception as e:
            keys = []

        for key in keys:
            label_dict[key] = self.h5.get('model/labels')[key].value

        plt.figure(figsize=(6,5))
        plt.plot(z,working_profile)
        valid = np.where(working_profile)[0]
        plt.xlim((z[valid[0]],z[valid[-1]]))
        plt.autoscale(False)
        for label in label_dict.keys():
            label_z = z[label_dict[label]]
            text_x = label_z
            text_y = working_profile[label_dict[label]]
            plt.text(text_x,text_y,label,ha='center',va='bottom')
        plt.xlabel('depth ($\mu$m)')
        plt.ylabel('linear intensity (ADU)')
        plt.title(tag)
        
        outfn = os.path.join(self.report_directory,tag.replace(' ','_')+'.png')
        plt.savefig(outfn)
        plt.close()


    def dispersion_report(self,N=None):
        if N is None:
            N = 20

        logs = []
        log_names = []

        try:
            for key in self.h5.get('dispersion/logs').keys():
                logs.append(self.h5.get('dispersion/logs/%s'%key))
                log_names.append(key)

        except:
            self.logger.info('Cannot retrieve dispersion logs from h5 store.')
            return

        try:
            coefs = self.h5.get('dispersion/coefficients')
        except Exception as e:
            self.logger.info('Cannot retrieve dispersion coefficients from h5 store.')

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
        plt.figure(figsize=(4*nlog,5))
        for idx,(log,log_name) in enumerate(zip(logs,log_names)):
            plt.subplot(1,nlog,1+idx)
            plot(log,log_name)

        outfn = os.path.join(self.report_directory,'dispersion_optimization.png')
        plt.savefig(outfn)
        plt.close()

    def processed_report(self,show=True):
        proc = np.abs(self.h5.get('processed_data')[:])
        n_volumes,n_slow,n_depth,n_fast = proc.shape
        for i_volume in range(n_volumes):
            outdir = os.path.join(self.report_directory,'processed_volume_%03d'%i_volume)
            self.makedirs(outdir)

            vol = proc[i_volume,:,:,:]

            test_vol = vol[5:-5:5,:,5:-5:5]
            cmin = np.median(test_vol)
            cmax = np.percentile(test_vol,99.95) # saturate 0.05% of pixels

            dpi = 100.0
            plt.figure(figsize=(n_fast/dpi,n_depth/dpi))
            
            for i_slow in range(n_slow):
                plt.cla()
                plt.axes([0,0,1,1])
                plt.imshow(vol[i_slow,:,:],cmap='gray',clim=(cmin,cmax),interpolation='none',aspect='auto')
                plt.savefig(os.path.join(outdir,'%04d.png'%i_slow),dpi=dpi)
                if show:
                    plt.pause(.0000001)

            plt.close()
        

if __name__=='__main__':

    h5 = H5('./oct_test_volume/oct_test_volume_2T.hdf5')
    r = Reporter(h5,'./oct_test_volume/oct_test_volume_2T_report')
    r.processed_report()
