import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import correlate2d,convolve2d


class StripRegistrar:

    def __init__(self,ref,strip_width=1):
        """Strips must be horizontal. If the strips to be registered are vertical,
        make sure to transpose before running this, and then transpose back when
        you're done."""
        self.ref = ref
        self.strip_width = strip_width

    def strip_register_image(self,im):
        sy,sx = im.shape

        y1s = np.arange(0,sy,self.strip_width).astype(np.int)
        y2s = y1s + self.strip_width
        y2s[-1] = min(y2s[-1],sy)


        xshifts = []
        yshifts = []
        corrs = []

        ref = (self.ref - np.mean(self.ref))/np.std(self.ref)
        
        for y1,y2 in zip(y1s,y2s):
            strip = im[y1:y2,:]
            strip = (strip - np.mean(strip))/np.std(strip)
            nxcval = np.real(np.fft.fftshift(np.fft.ifft2(np.fft.fft2(strip,s=(sy,sx))*np.conj(np.fft.fft2(self.ref)))))

            #nxcval = correlate2d(strip,ref)
            
            plt.cla()
            plt.imshow(nxcval)
            plt.pause(.1)

            peakVal = np.max(nxcval)
            peakIdx = np.where(nxcval==peakVal)

            yoff,xoff = peakIdx[0][0],peakIdx[1][0]
            
            if nxcval.shape[0]%2:
                yshift = (nxcval.shape[0]-1)/2.0 - yoff
            else:
                yshift = nxcval.shape[0]/2.0 - yoff

            yshift = yshift%sy

            if nxcval.shape[1]%2:
                xshift = (nxcval.shape[1]-1)/2.0 - xoff
            else:
                xshift = nxcval.shape[1]/2.0 - xoff

            xshifts.append(xshift)
            yshifts.append(yshift)
            corrs.append(peakVal/ref.shape[1]/)
            print corrs[-1]

        return y1s,y2s,yshifts,xshifts,corrs

if __name__=='__main__':


    template_height = 200
    template_width = 300
    im_height = 100
    im_width = 100
    
    template = np.zeros((template_height,template_width))

    for k in range(2000):
        x = np.random.randint(template_width)
        y = np.random.randint(template_height)
        template[y,x] = 1.0
    diameter = 10
    kernel = np.zeros((diameter,diameter))
    XX,YY = np.meshgrid(np.arange(diameter),np.arange(diameter))
    XX = XX - diameter/2.0
    YY = YY - diameter/2.0
    d = np.sqrt(XX**2+YY**2)
    kernel = np.zeros_like(d)
    kernel[np.where(d<=diameter/2.0)] = 1.0
    template = convolve2d(template,kernel,mode='same')

    
    ref = template[50:50+im_height,50:50+im_width]

    dx = 50
    dy = 50
    tar = np.zeros_like(ref)
    for k in range(100):
        dx = dx + np.random.randn()*0.0
        dy = dy + np.random.randn()*0.0
        tar[k,:] = template[k+round(dy),dx:dx+im_width]

    sr = StripRegistrar(ref,strip_width=5)
    print sr.strip_register_image(tar)

    
        
        
