import sys,os
import logging
import numpy as np
import scipy as sp
from scipy.ndimage.filters import generic_filter
from scipy.interpolate import bisplrep,bisplev
from matplotlib import pyplot as plt
from scipy.ndimage.morphology import grey_opening
from scipy.ndimage.filters import median_filter
from utils import translation,translation1,autotrim_bscan,find_peaks,shear,Clock,lateral_smooth_3d,polyfit2d,polyval2d
from octopod.Misc import H5
import octopod_config as ocfg
import logging
logging.basicConfig(level=logging.DEBUG)


class Cropper:

    def __init__(self,h5,debug=False):
        """Initialize model. May pass an h5py.File object or .hdf5 filename.
        The file or object must contain a 'processed_data' dataset containing
        at least one volume."""
        self.logger = logging.getLogger(__name__)
        if type(h5)==str:
            self.h5 = H5(h5)
            self.logger.info('Opening file %s.'%self.h5)
        else:
            self.h5 = h5



    def crop(self):

        vols = np.abs(self.h5.get('/processed_data')[:])

        nv,ns,nd,nf = vols.shape

        sumvol = np.zeros_like(vols[0,:,:,:])

        for k in range(nv):
            sumvol = sumvol + vols[k,:,:,:]


        self.profile = np.mean(np.mean(sumvol,axis=2),axis=0)

        global clicks
        crop_lims = self.autocrop(do_plot=True)
        clicks = [crop_lims[0],crop_lims[1]]
        print clicks
        
        fig = plt.figure(figsize=(10,6))
        
        def onclick(event):
            global clicks
            newclick = round(event.xdata)
            clicks.append(newclick)
            print clicks,
            if len(clicks)>=2:
                clicks = clicks[-2:]

            print clicks
            plt.cla()
            plt.plot(self.profile)
            for click in clicks:
                plt.axvline(click)
            plt.draw()
            
        cid = fig.canvas.mpl_connect('button_press_event',onclick)
        
        plt.plot(self.profile)
        if len(clicks):
            for click in clicks:
                plt.axvline(click)
                
        plt.show()

        if len(clicks)>=2:
            clicks = clicks[-2:]
            z1 = np.min(clicks)
            z2 = np.max(clicks)
            self.profile = self.profile[z1:z2+1]

        print self.h5.get('model').keys()


    def autocrop(self,noise_border=200,do_plot=False):

        x_in = np.arange(noise_border)
        y_in = self.profile[:noise_border]

        p = np.polyfit(x_in,np.log(y_in),deg=1)

        x_out = np.arange(len(self.profile))
        y_out_log = np.polyval(p,x_out)
        
        y_out = np.exp(y_out_log+.1)
        test = self.profile-y_out
        valid = np.where(test>0)[0]

        if do_plot:
            plt.figure()
            plt.plot(x_in,y_in,'ko',label='data for fit')
            plt.plot(x_out,self.profile,'g:',label='data')
            plt.plot(x_out,y_out,'r-',label='fit')
            plt.legend()
            plt.show()
        return valid[0],valid[-1]


def test():
    h5 = H5('./oct_test_volume/oct_test_volume_2T.hdf5')
    c = Cropper(h5,True)
    c.crop()
    
if __name__=='__main__':
    test()
