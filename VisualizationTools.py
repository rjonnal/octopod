from octopod.DataStore import H5
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import octopod_config as ocfg
import logging,sys


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
    #h5 = H5('./oct_test_volume/oct_test_volume_2T.hdf5')
    #pm = ProjectionMaker(h5)
    pm = ProjectionMaker('./oct_test_volume/oct_test_volume_2T.hdf5')
    
    pm.project()
