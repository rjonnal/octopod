import matplotlib
#matplotlib.use('Qt4Agg')
#matplotlib.interactive(True)
from octopod.DataStore import H5
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import octopod_config as ocfg
import logging,sys,os
#from mayavi import mlab
from utils import autotrim_volume, autotrim_bscan
import glob

# class Viewer3D:

#     def __init__(self,vol,zmin=None,zmax=None):
#         self.logger = logging.getLogger(__name__)
#         self.avol = vol
#         prof = np.mean(np.mean(self.avol,axis=2),axis=0)

#         if zmin is None or zmax is None:
#             fig = plt.figure()

#             global clicks,z1,z2
#             clicks = []
#             def onclick(event):
#                 global clicks,z1,z2
#                 newclick = round(event.xdata)
#                 clicks.append(newclick)
#                 if len(clicks)>=2:
#                     z1 = min(clicks[-2],clicks[-1])
#                     z2 = max(clicks[-2],clicks[-1])
#                     plt.cla()
#                     plt.plot(prof)
#                     plt.axvspan(z1,z2,alpha=0.3)
#                     plt.draw()

#             cid = fig.canvas.mpl_connect('button_press_event',onclick)
#             plt.plot(prof)
#             plt.show()
#             zmin = z1
#             zmax = z2

#         self.avol = self.avol[:,zmin:zmax,:]
        
#         vmed = np.median(self.avol)
#         vmean = np.mean(self.avol)
#         vmax = np.max(self.avol)
#         vmin = np.min(self.avol)
#         vstd = np.std(self.avol)
        
#         scatter_field = mlab.pipeline.scalar_field(self.avol)
#         vmin,vmax = np.percentile(self.avol,[25,95])
#         mlab.pipeline.volume(scatter_field,vmin=vmin,vmax=vmax,color='gray')
#         # mlab.pipeline.iso_surface(scatter_field, contours=[vmax, ],)
#         # mlab.pipeline.image_plane_widget(scatter_field,
#         #                     plane_orientation='z_axes',
#         #                     slice_index=10)
        
#         mlab.show()

class VolumeProjectionMaker:

    def __init__(self,avol,outdir='./vpm_saves'):
        self.logger = logging.getLogger(__name__)

        self.outdir = outdir
        self.logger.info('Initializing VolumeProjectionMaker.')
        self.avol = avol
        self.logger.info('Computing slow-axis projection.')
        self.slow_projection = np.mean(self.avol,axis=0)
        self.logger.info('Computing fast-axis projection.')
        self.fast_projection = np.mean(self.avol,axis=2).T
        self.logger.info('Computing depth profile.')
        self.profile = np.mean(self.fast_projection,axis=1)
        outdir_files = glob.glob(os.path.join(self.outdir,'*.*'))
        self.logger.info('Existing files:')
        for f in outdir_files:
            self.logger.info(f)

        
    def project(self):

        global clicks
        clicks = []
        
        fig = plt.figure(figsize=(10,6))
        
        def onpress(event):
            if len(clicks)>=2:
                if event.key=='enter':
                    if not os.path.exists(self.outdir):
                        self.logger.info('Making %s.'%self.outdir)
                        os.makedirs(self.outdir)
                    z1 = min(clicks[-2],clicks[-1])
                    z2 = max(clicks[-2],clicks[-1])
                    label = raw_input('Please type a label for this projection: ')
                    tag = '%s_%d_%d'%(label,z1,z2)
                    png_fn = os.path.join(self.outdir,'%s.png'%tag)
                    #npy_fn = os.path.join(self.outdir,'%s.npy'%tag)
                    self.logger.info('Saving PNG to %s.'%png_fn)
                    plt.savefig(png_fn)

        cid = fig.canvas.mpl_connect('key_press_event',onpress)
                    
        def onclick(event):
            global clicks
            newclick = round(event.xdata)
            clicks.append(newclick)
            if len(clicks)>=2:
                z1 = int(np.round(min(clicks[-2],clicks[-1])))
                z2 = int(np.round(max(clicks[-2],clicks[-1])))
                self.areal_projection = np.mean(self.avol[:,z1:z2,:],axis=1)
                self.clim = np.percentile(self.areal_projection[np.where(self.areal_projection)],(1,99.5))
                plt.subplot(1,2,1)
                plt.cla()
                plt.plot(self.profile)
                plt.axvspan(z1,z2,alpha=0.3)
                plt.subplot(1,2,2)
                plt.cla()
                #plt.imshow(self.areal_projection,aspect='auto',interpolation='none')
                plt.imshow(self.areal_projection,aspect='auto',interpolation='none',cmap='gray',clim=self.clim)
                plt.draw()
            
        cid = fig.canvas.mpl_connect('button_press_event',onclick)


        self.areal_projection = np.mean(self.avol,axis=1)
        self.clim = np.percentile(self.areal_projection[np.where(self.areal_projection)],(1,99.5))
        plt.subplot(1,2,1)
        plt.plot(self.profile)
        plt.subplot(1,2,2)
        plt.imshow(self.areal_projection,aspect='auto',interpolation='none',cmap='gray',clim=self.clim)
        plt.show()





if __name__=='__main__':
    fn = './oct_test_volume/oct_test_volume_2T.hdf5'
    #pm = ProjectionMaker(fn)
    #pm.project()

    v3d = Viewer3D(fn)
