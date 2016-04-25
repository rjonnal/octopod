import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import correlate2d,convolve2d,fftconvolve
from scipy.misc import imread
import sys
from octopod import utils
from time import sleep
import logging
logging.basicConfig(level='INFO')


class StripRegistrar:
    def __init__(self,ref,strip_width=1,do_plots=False):
        """Strips must be horizontal. If the strips to be registered are vertical,
        make sure to transpose before running this, and then transpose back when
        you're done."""
        self.ref = ref
        self.strip_width = strip_width
        self.rsy,self.rsx = self.ref.shape
        self.do_plots = do_plots
        self.strips = []
        self.y1s = []
        self.y2s = []
        self.x1s = []
        self.x2s = []
        self.layers = []
        self.corrs = []
        self.current_layer = 0
        self.logger = logging.getLogger(__name__)
        self.logger.info('Creating StripRegistrar object.')
        
    def add(self,im,xmax=None,ymax=None):
        y1s,y2s,yshifts,xshifts,corrs = self.get_reg_info(im,xmax=xmax,ymax=ymax)
        self.strips = self.strips + [im[y1:y2,:] for y1,y2 in zip(y1s,y2s)]
        self.y1s = self.y1s + yshifts
        self.y2s = self.y2s + [ytemp + self.strip_width for ytemp in yshifts]
        self.x1s = self.x1s + xshifts
        self.x2s = self.x2s + [xtemp + im.shape[1] for xtemp in xshifts]
        self.layers = self.layers + [self.current_layer]*len(y1s)
        self.corrs = self.corrs + corrs
        self.current_layer = self.current_layer + 1

    def render(self,corr_percentile=50):
        corr_thresh = np.percentile(self.corrs,corr_percentile)
        self.logger.info('render: correlation threshold: %0.2f'%corr_thresh)
        valid_idx = np.where(self.corrs>=corr_thresh)

        def valid(lst):
            return [item for item,corr in zip(lst,self.corrs) if corr>=corr_thresh]
        
        y1s = valid(self.y1s)
        y2s = valid(self.y2s)
        x1s = valid(self.x1s)
        x2s = valid(self.x2s)
        strips = valid(self.strips)

        yoffset = -np.min(y1s)
        xoffset = -np.min(x1s)

        y1s = [y1-np.min(y1s) for y1 in y1s]
        y2s = [y1 + s.shape[0] for y1,s in zip(y1s,strips)]

        x1s = [x1-np.min(x1s) for x1 in x1s]
        x2s = [x1 + s.shape[1] for x1,s in zip(x1s,strips)]

        if yoffset>=0:
            # all strips are below top of reference
            ry1 = 0
            ry2 = self.rsy
            outsy = max(np.max(y2s)+yoffset,ry2)
        if xoffset>=0:
            # all strips are to the right of reference
            rx1 = 0
            rx2 = self.rsx
            outsx = max(np.max(x2s)+xoffset,rx2)
        if yoffset<0:
            # some strips are above top of reference
            ry1 = -yoffset
            ry2 = ry1 + self.rsy
            outsy = max(np.max(y2s)+yoffset,ry2)
        if xoffset<0:
            # some strips start to the left of reference left edge
            rx1 = -xoffset
            rx2 = rx1 + self.rsx
            outsx = max(np.max(x2s)+xoffset,rx2)

        sumimage = np.zeros((outsy,outsx))
        #sumimage[ry1:ry2,rx1:rx2] = self.ref
        #utils.scaleshow(sumimage)
        counterimage = np.zeros_like(sumimage)
        print xoffset,yoffset,ry1,ry2
        
        y1s = [y1 - yoffset for y1 in y1s]
        y2s = [y2 - yoffset for y2 in y2s]
        x1s = [x1 + xoffset for x1 in x1s]
        x2s = [x2 + xoffset for x2 in x2s]
        
        #sumimage = np.zeros((np.max(y2s),np.max(x2s)))
        #counterimage = np.ones_like(sumimage)

        for y1,y2,x1,x2,strip in zip(y1s,y2s,x1s,x2s,strips):
            while(y1<0):
                self.logger.info('render: cropping strip')
                strip = strip[1:,:]
                y1 = y1 + 1

            while (y2>outsy):
                self.logger.info('render: cropping strip')
                strip = strip[:-1,:]
                y2 = y2 - 1

            while (x1<0):
                self.logger.info('render: cropping strip')
                strip = strip[:,1:]
                x1 = x1 + 1

            while (x2>outsx):
                self.logger.info('render: cropping strip')
                strip = strip[:,:-1]
                x2 = x2 - 1
            
            infostring = 'inserting %d x %d strip at %d, %d'%(strip.shape[0],strip.shape[1],y1,x1)
            self.logger.info('render: %s'%infostring)
            sumimage[y1:y2,x1:x2] = sumimage[y1:y2,x1:x2] + strip
            counterimage[y1:y2,x1:x2] = counterimage[y1:y2,x1:x2] + 1

        #utils.scaleshow(counterimage)
        #utils.scaleshow(sumimage/counterimage,clim=utils.scaleshow(self.ref))
        counterimage[np.where(counterimage==0)]=1.0e-9
        return sumimage,counterimage
        
        
        

    def dewarp_image(self,im,reg_info,corr_threshold=0.0):
        y1s = np.array(reg_info[0])
        y2s = np.array(reg_info[1])
        yshifts = np.array(reg_info[2])
        xshifts = np.array(reg_info[3])
        corrs = np.array(reg_info[4])

        valid_idx = np.where(corrs>=corr_threshold)[0]

        y1s = y1s[valid_idx]
        y2s = y2s[valid_idx]
        yshifts = yshifts[valid_idx]
        xshifts = xshifts[valid_idx]
        
        ysmin = np.min(yshifts)
        xsmin = np.min(xshifts)
        
        yshifts = yshifts - ysmin
        xshifts = xshifts - xsmin
        plt.plot(yshifts)
        plt.show()
        sy,sx = im.shape
        out = np.zeros((np.max(yshifts)+self.strip_width,sx+np.max(xshifts)))
        for y1,y2,ys,xs,corr in zip(y1s,y2s,yshifts,xshifts,corrs):
            print y1
            im_in = im[y1:y2,:]
            oy1 = ys
            oy2 = ys+self.strip_width
            oy2 = min(oy2,out.shape[0])
            ox1 = xs
            ox2 = xs+sx
            out[oy1:oy2,ox1:ox2] = im_in
            if self.do_plots:
                plt.figure(1)
                plt.cla()
                plt.imshow(im_in,interpolation='none')
                plt.figure(2)
                plt.cla()
                plt.imshow(out,interpolation='none')
                plt.pause(.1)
        return out

        
        
    def get_reg_info(self,im,xmax=None,ymax=None):

        # check to see if this is the reference image:
        epsilon = 1e-9
        test = np.abs(np.sum(im-self.ref))
        print test
        is_reference = test<epsilon
        if is_reference:
            self.logger.info('get_reg_info: this is the reference image')
        
        sy,sx = im.shape

        y1s = np.arange(0,sy,self.strip_width).astype(np.int)
        y2s = y1s + self.strip_width
        y2s[-1] = min(y2s[-1],sy)

        xshifts = []
        yshifts = []
        corrs = []

        ref = (self.ref - np.mean(self.ref))/np.std(self.ref)

        if self.do_plots:
            plt.figure(figsize=(8,12))
        
        for idx,(y1,y2) in enumerate(zip(y1s,y2s)):
            self.logger.info('get_reg_info: registering strip %d of %d.'%(idx+1,len(y1s)))
            expected_y = self.rsy - y1
            expected_x = self.rsx

            if is_reference:
                yshift = self.rsy - expected_y - 2
                xshift = self.rsx - expected_x - 2
                if xshift>sx/2:
                    xshift = xshift - sx + 2
                xshifts.append(xshift)
                yshifts.append(yshift)
                corrs.append(-np.inf)
                continue
                
            
            strip = im[y1:y2,:]
            strip = (strip - np.mean(strip))/np.std(strip)
            nxcval = fftconvolve(strip,self.ref[::-1,::-1],mode='full')
            #nxcval = np.real(np.fft.fftshift(np.fft.ifft2(np.fft.fft2(strip,s=(2*sy,2*sx))*np.conj(np.fft.fft2(self.ref,s=(2*sy,2*sx))))))
            #nxcval = np.real((np.fft.ifft2(np.fft.fft2(strip,s=(sy,sx))*np.conj(np.fft.fft2(self.ref,s=(sy,sx))))))

            if ymax is not None:
                ycrop1 = expected_y - ymax
                ycrop2 = expected_y + ymax
                ycrop1 = max(ycrop1,0)
                ycrop2 = min(ycrop2,nxcval.shape[0])
                newnxcval = np.zeros_like(nxcval)
                newnxcval[ycrop1:ycrop2,:] = nxcval[ycrop1:ycrop2,:]
                nxcval = newnxcval
                del newnxcval
                
            if xmax is not None:
                xcrop1 = expected_x - xmax
                xcrop2 = expected_x + xmax
                xcrop1 = max(xcrop1,0)
                xcrop2 = min(xcrop2,nxcval.shape[1])
                newnxcval = np.zeros_like(nxcval)
                newnxcval[:,xcrop1:xcrop2] = nxcval[:,xcrop1:xcrop2]
                nxcval = newnxcval
                del newnxcval

            if self.do_plots:
                plt.subplot(3,1,1)
                plt.cla()
                plt.imshow(im,interpolation='none',cmap='gray')
                plt.autoscale(False)
                plt.axhspan(y1,y2,alpha=0.3)
                plt.subplot(3,1,2)
                plt.cla()
                plt.imshow(self.ref,interpolation='none',cmap='gray')
                plt.subplot(3,1,3)
                plt.cla()
                plt.imshow(nxcval)
                plt.autoscale(False)
                py,px = np.where(nxcval==np.max(nxcval))
                plt.plot(px,py,'ks')
                plt.pause(.1)

            peak_object = np.where(nxcval==np.max(nxcval))

            ypeak = peak_object[0][0]
            xpeak = peak_object[1][0]
            
            yshift = self.rsy - ypeak - 2
            xshift = self.rsx - xpeak - 2
            if xshift>sx/2:
                xshift = xshift - sx + 2
            corr = np.max(nxcval)/self.rsx/self.rsy/self.strip_width
            xshifts.append(xshift)
            yshifts.append(yshift)
            corrs.append(corr)
            
        if self.do_plots:
            plt.close()


        cmax = np.max(corrs)
        newcorrs = []
        for corr in corrs:
            if corr==-np.inf:
                newcorrs.append(cmax)
            else:
                newcorrs.append(corr)
        corrs = newcorrs

        # for y1,y2,ys,xs,corr in zip(y1s,y2s,yshifts,xshifts,corrs):
        #     print y1,y2,ys,xs,corr

        # sys.exit()
        
        return y1s,y2s,yshifts,xshifts,corrs

if __name__=='__main__':
    
    template_height = 400
    template_width = 800
    im_height = 200
    im_width = 600
    strip_width = 4
    xoff = 50
    yoff = 50
    
    template = imread('./invivo_b.png')
    if len(template.shape)==3:
        template = np.mean(template,axis=2)


    ref = template[yoff:yoff+im_height,xoff:xoff+im_width]
    def make_tar(mode='random'):
        dx = xoff+.5
        dy = yoff+.5
        tar = np.zeros_like(ref)
        for k in range(im_height):
            rdx = round(dx)
            rdy = round(dy)
            tar[k,:] = template[k+rdy,rdx:rdx+im_width]
            if mode=='linear_shear':
                dx = dx + 1.0/strip_width
                dy = dy
            elif mode=='random':
                dx = dx + np.random.randn()*1.0
                dy = dy + np.random.randn()*1.0
            elif mode=='no_shift':
                pass
        return tar

    sr = StripRegistrar(ref,strip_width=strip_width,do_plots=False)
    
    for k in range(10):
        print k
        sr.add(make_tar(),xmax=20,ymax=50)
        
    registered = sr.render(90)
    plt.show()
    

    
        
        
