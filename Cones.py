import numpy as np
from scipy import ndimage
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level='INFO')

class Mosaic:

    def __init__(self,im):
        self.cones = im
        self.logger = logging.getLogger(__name__)
        self.logger.info('Creating Mosaic object.')
        
    def get_coordinates(self,left_crop=0,show_plot=True,minimum_distance=0,png_filename=None,clim=None):
        self.logger.info('get_coordinates: starting.')
        cones = self.cones
        cones[np.where(np.isnan(cones))] = 0

        coness = ndimage.gaussian_filter(cones,1.1)
        conesf = ndimage.gaussian_filter(cones,.75)

#        conesf = ndimage.gaussian_filter(cones,0.5)

        self.logger.info('get_coordinates: computing derivatives.')
        conesdv = np.diff(conesf,axis=0)
        top = conesdv[:-1,:]
        bottom = conesdv[1:,:]
        tbDirectionMask = np.zeros(top.shape)
        tbDirectionMask[np.where(top>bottom)] = 1
        

        conesdh = np.diff(conesf,axis=1)
        left = conesdh[:,:-1]
        right = conesdh[:,1:]
        lrDirectionMask = np.zeros(left.shape)
        lrDirectionMask[np.where(left>right)] = 1
        
        hmap = (left*right*(-1)).clip(0,1)*lrDirectionMask
        vmap = (top*bottom*(-1)).clip(0,1)*tbDirectionMask
        hmapfull = np.zeros(cones.shape)
        vmapfull = np.zeros(cones.shape)
        hmapfull[:,1:-1] = hmap
        vmapfull[1:-1,:] = vmap
        hmap = hmapfull
        vmap = vmapfull
        cmap = hmap*vmap
        cmap[:,:left_crop] = 0

        lowIntensityMask = np.ones(conesf.shape)
        thresh = conesf.mean()# + conesf.std()
        lowIntensityMask[np.where(conesf<thresh)] = 0
#        lowIntensityMask = lowIntensityMask[1:-1,1:-1]
        cmap = cmap * lowIntensityMask


        cy,cx = cmap.nonzero()


        if minimum_distance:
            self.logger.info('get_coordinates: checking distances.')
            nCones = len(cy)
            valid = np.ones(nCones)
            for i in range(nCones):
                for j in range(i+1,nCones):
                    d = np.sqrt((cy[i]-cy[j])**2+(cx[i]-cx[j])**2)
                    if d<minimum_distance:
                        valid[j] = 0
            cy = cy[np.where(valid)]
            cx = cx[np.where(valid)]
                

        if show_plot or png_filename is not None:

            if clim is None:
                clim = np.percentile(conesf,[5,99.5])
            
            self.logger.info('get_coordinates: making plot.')
            plt.figure(figsize=(12,6),facecolor='black')
            #plt.subplot(2,4,1)
            plt.axes((0,.5,.2475,.485),axisbg='black')
            imh = plt.imshow(conesf,aspect='auto')

            imh.set_clim(clim)

            plt.xticks([])
            plt.yticks([])
            #plt.title('original image')
            #plt.subplot(2,4,2)
            plt.axes((.25,.5,.2475,.485),axisbg='black')
            imh = plt.imshow(coness,aspect='auto')
            imh.set_clim(clim)


            plt.xticks([])
            plt.yticks([])
            #plt.title('gaussian convolved image')
            #plt.subplot(2,4,3)
            plt.axes((.5,.5,.2475,.485),axisbg='black')
            plt.imshow(conesdv,aspect='auto')
            plt.xticks([])
            plt.yticks([])
            #plt.title('vertical gradient')
            #plt.subplot(2,4,4)
            plt.axes((.75,.5,.2475,.485),axisbg='black')
            plt.imshow(conesdh,aspect='auto')
            plt.xticks([])
            plt.yticks([])
            #plt.title('horizontal gradient')
            #plt.subplot(2,4,5)
            plt.axes((0,0,.2475,.485),axisbg='black')
            plt.imshow(vmap,aspect='auto')
            plt.xticks([])
            plt.yticks([])
            #plt.title('vertical falling zeros')
            #plt.subplot(2,4,6)
            plt.axes((.25,0,.2475,.485),axisbg='black')
            plt.imshow(hmap,aspect='auto')
            plt.xticks([])
            plt.yticks([])
            #plt.title('horizontal falling zeros')
            #plt.subplot(2,4,7)
            plt.axes((.5,0,.2475,.485),axisbg='black')
            plt.imshow(cmap,aspect='auto')
            plt.set_cmap('gray')
            plt.xticks([])
            plt.yticks([])
            #plt.title('intersection of ridges')
            #plt.subplot(2,4,8)
            plt.axes((0.75,0,.2475,.485),axisbg='black')

            imh = plt.imshow(conesf)
            imh.set_clim(clim)
            plt.xticks([])
            plt.yticks([])
            ph = plt.plot(cx,cy,'go')[0]
            ph.set_markersize(4)
            
            plt.axis('tight')

            if png_filename is not None:
                plt.savefig(png_filename,dpi=300,facecolor='k',edgecolor='none')

            if show_plot:
                plt.show()
            else:
                plt.close()
                
        return cy,cx
