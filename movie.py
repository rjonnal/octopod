import numpy as np
from matplotlib import pyplot as plt
import scipy as sp
import png
import datetime
import os
import sys
from subprocess import call
import shutil
import matplotlib as mpl
import platform
from glob import glob

class Movie:

    def __init__(self,avifn=None,fps=10,cmin=None,cmax=None,dpi=96,make_avi=True,make_webm=True,make_wmv=True,autoclean=True,thumbfn=None,previewfn=None):
        self._osname = platform.system().lower()
        self._fps = fps
        self._cmin = cmin
        self._cmax = cmax
        self._dpi = dpi
        self._make_webm = make_webm
        self._make_wmv = make_wmv
        self._make_avi = make_avi
        self._autoclean = autoclean
        
        nowStr = datetime.datetime.now().strftime("%Y%m%d%H%M%S")

        if avifn is None:
            self._avifn = os.path.join('.','%s.avi'%nowStr)
        else:
            self._avifn = avifn

        if thumbfn is not None:
            self._thumbfn = thumbfn
        else:
            self._thumbfn = self._avifn.replace('.avi','_thumb.png')
        
        if previewfn is not None:
            self._previewfn = previewfn
        else:
            self._previewfn = self._avifn.replace('.avi','_preview.png')

        self._webmfn = self._avifn.replace('.avi','.webm')
        self._wmvfn = self._avifn.replace('.avi','.wmv')

        #self._wdir = './tmp_%s%0.5f'%(nowStr,np.random.rand())
        cwd = os.getcwd()
        self._wdir = os.path.join(cwd,'tmp_%s%0.5f'%(nowStr,np.random.rand()))
        os.makedirs(self._wdir)

        self._counter = 0

    def __del__(self):
        if os.path.exists(self._wdir) and self._autoclean:
            try:
                print 'Movie: abnormal exit, deleting %s.'%self._wdir
                shutil.rmtree(self._wdir)
            except Exception as e:
                print e

    def add(self,im):
        plt.ion()
        fulloutfn = os.path.join(self._wdir,'frame_%09d.png'%self._counter)
        if type(im)==np.ndarray:
            if self._cmin is None:
                cmin = im.min()
            else:
                cmin = self._cmin

            if self._cmax is None:
                cmax = im.max()
            else:
                cmax = self._cmax

            bmp = np.uint8(np.round((im-cmin)/(cmax-cmin)*255).clip(0,255))
            bmin = bmp.min()
            bmax = bmp.max()
            png.from_array(bmp,'L').save(fulloutfn)

            if self._counter<5:
                thumb = sp.ndimage.zoom(bmp,.2)
                png.from_array(thumb,'L').save(self._thumbfn)
                png.from_array(bmp,'L').save(self._previewfn)

            if self._counter==0:
                plt.figure()
                self._imh = plt.imshow(bmp,interpolation='none')
                plt.set_cmap('gray')
                self._imh.set_clim((bmin,bmax))
                plt.draw()
                
            else:
                self._imh.set_data(bmp)
                self._imh.set_clim((bmin,bmax))

            plt.draw()

        elif type(im)==mpl.figure.Figure:
            im.savefig(fulloutfn,dpi=self._dpi)
            if self._counter==0:
                plt.show()
            else:
                plt.draw()

        self._counter = self._counter + 1
        plt.ioff()

    def make(self):
        plt.close()
        if self._osname == 'linux':
            if self._make_avi:
                command = ['mencoder', 'mf://%s'%os.path.join(self._wdir,'frame*.png'), '-mf', 'type=png:fps=%d'%self._fps,
                        '-ovc','lavc','-o',self._avifn]
                print command
                call(command)

            if self._make_webm:

                command = ['mencoder', 'mf://%s'%os.path.join(self._wdir,'frame*.png'), '-mf', 'type=png:fps=%d'%self._fps,
                          '-ovc','lavc', '-lavcopts', 'threads=4:acodec=vorbis:vcodec=libvpx', '-ffourcc', 'VP80','-o',self._webmfn]
                call(command)

            if self._make_wmv:

                command = ['mencoder', 'mf://%s'%os.path.join(self._wdir,'frame*.png'), '-mf', 'type=png:fps=%d'%self._fps,
                           '-ovc','lavc', '-lavcopts', 'threads=4:vbitrate=10000:acodec=wmav1:vcodec=wmv1', '-o',self._wmvfn]
                call(command)

        if self._osname == 'windows':
            if self._make_webm:
                fstring = os.path.join(self._wdir,'frame_%09d.png')
                command = ['ffmpeg', '-y', '-r', '%d'%(self._fps), '-i', '%s'%fstring,
                                '-r', '%d'%(self._fps), '-crf', '4', '-b:v', '256k', self._webmfn]
                call(command,shell=True)
                
        if self._autoclean:
            print 'Movie: normal exit, trying to delete %s... '%self._wdir,
            try:
                shutil.rmtree(self._wdir)
                print 'success!'
            except Exception as e:
                print 'uh oh: %s'%e
            
