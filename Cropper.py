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
        """Initialize Cropper. May pass an h5py.File object or .hdf5 filename.
        The cropper crops the processed volume data cube, interactively, and
        updates all of the following, when they exist, because each of them
        is defined (implicitly) relative to the coordinates of the original
        processed data cube:
        model profile
        model labels
        z offsets and z offset fit
        volume labels
        """
        self.logger = logging.getLogger(__name__)
        if type(h5)==str:
            self.h5 = H5(h5)
            self.logger.info('Opening file %s.'%self.h5)
        else:
            self.h5 = h5

    def crop(self):
        self.logger.info('crop: Starting')
        vols = np.abs(self.h5.get('/processed_data')[:])

        nv,ns,nd,nf = vols.shape

        sumvol = np.zeros_like(vols[0,:,:,:])

        for k in range(nv):
            sumvol = sumvol + vols[k,:,:,:]


        self.profile = np.mean(np.mean(sumvol,axis=2),axis=0)

        global clicks

        if len(self.profile)==1024:
            crop_lims = self.autocrop(do_plot=False)
            clicks = [crop_lims[0],crop_lims[1]]
        else:
            clicks = []
            
        fig = plt.figure(figsize=(10,6))
        
        def onclick(event):
            global clicks
            newclick = round(event.xdata)
            clicks.append(newclick)
            if len(clicks)>=2:
                clicks = clicks[-2:]

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


        plt.title('Click crop points z1, z2 in order and close window to complete cropping.')
        plt.show()

        if len(clicks)==2:
            clicks = clicks[-2:]
        else:
            self.logger.info('crop: Oops %d points selected. Exiting.'%len(clicks))
            return

        z1 = clicks[0]
        z2 = clicks[1]

        self.profile = self.profile[z1:z2]

        def write(loc,val):
            self.logger.info('crop: Writing to %s'%loc)
            try:
                self.h5.put(loc,val)
            except Exception as e:
                self.logger.error('crop: Error: could not complete write.')

        try:
            # crop the model profile:
            model_profile = self.h5['/model/profile'][:]
            model_profile = model_profile[z1:z2]
            write('/model/profile',model_profile)
        except Exception as e:
            self.logger.info('Cannot crop model profile: %s'%e)

        try:
            # fix the model labels:
            for key in self.h5['/model/labels/'].keys():
                val = self.h5['/model']['labels'][key].value
                val = val - z1
                write('/model/labels/%s'%key,val)
        except Exception as e:
            self.logger.info('Cannot fix model labels: %s'%e)

        try:
            # fix the 2D layer labels
            for key in self.h5['/model/volume_labels'].keys():
                val = self.h5['/model']['volume_labels'][key][:]
                val = val - z1
                write('/model/volume_labels/%s'%key,val)
        except Exception as e:
            self.logger.info('Cannot fix volume labels: %s'%e)

        try:
            # fix z-offsets (offset between model and processed data cube)
            val = self.h5['/model/z_offsets'][:]
            val = val - z1
            write('/model/z_offsets',val)
        except Exception as e:
            self.logger.info('Cannot fix offsets: %s'%e)

        try:
            # fix the offset fit (a smoothed version of z-offsets)
            val = self.h5['/model/z_offset_fit'][:]
            val = val - z1
            write('/model/z_offset_fit',val)
        except Exception as e:
            self.logger.info('Cannot fix offset fit: %s'%e)

        try:
            # last, but not least, crop the data
            val = self.h5['/processed_data'][:]
            val = val[:,:,z1:z2,:]
            plt.figure()
            plt.imshow(np.abs(val[0,0,:,:]))
            plt.title('Ctrl-C to stop')
            plt.show()
            write('/processed_data',val)
        except Exception as e:
            self.logger.info('Cannot crop processed_data: %s'%e)

        self.h5.repack()
            
    def autocrop(self,noise_border=200,do_plot=False):

        x_in = np.arange(noise_border)
        y_in = self.profile[:noise_border]


        noise_rms = np.std(y_in)
        
        p = np.polyfit(x_in,np.log(y_in),deg=1)

        x_out = np.arange(len(self.profile))
        y_out_log = np.polyval(p,x_out)
        
        y_out = np.exp(y_out_log)+3*noise_rms


        
        test = self.profile-y_out

        
        valid = np.where(test>0)[0]

        if do_plot:
            plt.figure()
            plt.plot(x_in,y_in,'ko',label='data for fit')
            plt.plot(x_out,self.profile,'g:',label='data')
            plt.plot(x_out,y_out,'r-',label='fit')
            plt.legend()
            plt.show()
        return valid[0]-10,valid[-1]+10


def test():
    h5 = H5('./oct_test_volume/oct_test_volume_2T.hdf5')
    #m = Model(h5)
    #m.click_label()
    c = Cropper(h5,True)
    c.crop()
    
if __name__=='__main__':
    test()
