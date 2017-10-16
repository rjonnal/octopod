import sys
import matplotlib
from octopod.DataStore import H5
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import octopod_config as ocfg
import logging,sys,os
from mayavi import mlab
#mlab.options.backend = 'envisage'

from utils import autotrim_volume, autotrim_bscan
import glob
from tvtk.util.ctf import ColorTransferFunction
from tvtk.util.ctf import PiecewiseFunction

class Viewer3D:

    def __init__(self,vol,figsize=(800,600)):
        self.logger = logging.getLogger(__name__)
        self.data = vol
        
        self.vmed = np.median(self.data)
        self.vmean = np.mean(self.data)
        self.vmax = np.max(self.data)
        self.vmin = np.min(self.data)
        self.vstd = np.std(self.data)
        
        self.field = mlab.pipeline.scalar_field(self.data)
        self.volume = mlab.pipeline.volume(self.field)
        mlab.close()
        self.fig = mlab.figure(size=figsize)
        self.size =figsize
        
    def savefig(self,fn,size=None,magnification='auto'):
        if size is None:
            size = self.size
        mlab.savefig(fn,size=size,magnification=magnification,figure=self.fig)

    def screenshot(self,figure=None,mode='rgb',antialiased=False):
        return mlab.screenshot(figure,mode,antialiased)
            
    def show(self):
        mlab.show()
        sys.exit()

    def slice(self,plane_orientation,index):
        s = mlab.pipeline.image_plane_widget(self.field,plane_orientation=plane_orientation,slice_index=index)
        return s

    def set_colormap(self):
        values = np.linspace(self.vmin,self.vmax,256)
        lums = np.linspace(0,1,256)
        ctf = ColorTransferFunction()
        for v,l in zip(values,lums):
            ctf.add_rgb_point(v,l,l,l)
        self.volume._volume_property.set_color(ctf)
        self.volume._ctf = ctf
        self.volume.update_ctf = True



        # # Changing the ctf:
        # ctf.add_rgb_point(value, r, g, b)  # r, g, and b are float
        #                                    # between 0 and 1
        # ctf.add_hsv_point(value, h, s, v)
        # # ...
        # vol._volume_property.set_color(ctf)
        # vol._ctf = ctf
        # vol.update_ctf = True
    
        # # Changing the otf:
        # otf = PiecewiseFunction()
        # otf.add_point(value, opacity)
        # vol._otf = otf
        # vol._volume_property.set_scalar_opacity(otf)


        # mlab.pipeline.iso_surface(field, contours=[vmax, ],)
        # mlab.pipeline.image_plane_widget(scatter_field,
        #                     plane_orientation='z_axes',
        #                     slice_index=10)
        
        mlab.show()

