import matplotlib
matplotlib.use('Qt4Agg')
#matplotlib.interactive(True)

from octopod.DataStore import H5
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import octopod_config as ocfg
import logging,sys
from mayavi import mlab
from utils import autotrim_volume, autotrim_bscan

class Viewer3D:

    def __init__(self,h5,vidx=0,depth=128):
        self.logger = logging.getLogger(__name__)
        if type(h5)==str:
            h5 = H5(h5)

        self.avol = np.abs(h5.get('/processed_data'))[vidx,:,:,:]
        prof = np.mean(np.mean(self.avol,axis=2),axis=0)

        fig = plt.figure()
        
        global clicks,z1,z2
        clicks = []
        def onclick(event):
            global clicks,z1,z2
            newclick = round(event.xdata)
            clicks.append(newclick)
            if len(clicks)>=2:
                z1 = min(clicks[-2],clicks[-1])
                z2 = max(clicks[-2],clicks[-1])
                plt.cla()
                plt.plot(prof)
                plt.axvspan(z1,z2,alpha=0.3)
                plt.draw()
            
        cid = fig.canvas.mpl_connect('button_press_event',onclick)
        plt.plot(prof)
        plt.show()


        self.avol = self.avol[:,z1:z2,:]
        #self.avol = autotrim_volume(self.avol,depth)



        
        vmed = np.median(self.avol)
        vmean = np.mean(self.avol)
        vmax = np.max(self.avol)
        vmin = np.min(self.avol)
        vstd = np.std(self.avol)
        
        
        # plt.hist(np.ravel(self.avol),bins=100)
        # plt.axvline(vmean+vstd,color='r')
        # plt.axvline(vmean+2*vstd,color='r')
        # plt.axvline(vmean+3*vstd,color='r')
        # plt.axvline(vmean,color='g')
        # plt.show()
        scatter_field = mlab.pipeline.scalar_field(self.avol)
        mlab.pipeline.volume(scatter_field,vmin=vmean,vmax=vmed+3*vstd)
        mlab.show()

class ProjectionMaker:

    def __init__(self,h5,vidx=0):
        self.logger = logging.getLogger(__name__)
        if type(h5)==str:
            h5 = H5(h5)
        self.logger.info('Extracting volume from hdf5 file.')
        self.avol = np.abs(h5.get('/processed_data'))[vidx,:,:,:]
        self.logger.info('Computing slow-axis projection.')
        self.slow_projection = np.mean(self.avol,axis=0)
        self.logger.info('Computing fast-axis projection.')
        self.fast_projection = np.mean(self.avol,axis=2).T
        self.logger.info('Computing depth profile.')
        self.profile = np.mean(self.fast_projection,axis=1)
        
    def project(self):

        global clicks
        clicks = []
        
        fig = plt.figure(figsize=(10,6))
        
        def onclick(event):
            global clicks
            newclick = round(event.xdata)
            clicks.append(newclick)
            if len(clicks)>=2:
                z1 = min(clicks[-2],clicks[-1])
                z2 = max(clicks[-2],clicks[-1])
                self.areal_projection = np.mean(self.avol[:,z1:z2,:],axis=1)
                plt.subplot(1,2,1)
                plt.cla()
                plt.plot(self.profile)
                plt.axvspan(z1,z2,alpha=0.3)
                plt.subplot(1,2,2)
                plt.cla()
                plt.imshow(self.areal_projection,aspect='auto',interpolation='none')
                plt.draw()
            
        cid = fig.canvas.mpl_connect('button_press_event',onclick)


        self.areal_projection = np.mean(self.avol,axis=1)
        plt.subplot(1,2,1)
        plt.plot(self.profile)
        plt.subplot(1,2,2)
        plt.imshow(self.areal_projection,aspect='auto',interpolation='none')
        plt.show()
        





if __name__=='__main__':
    fn = './oct_test_volume/oct_test_volume_2T.hdf5'
    #pm = ProjectionMaker(fn)
    #pm.project()

    v3d = Viewer3D(fn)
