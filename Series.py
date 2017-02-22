import numpy as np
from matplotlib import pyplot as plt
import sys,os
from octopod import H5,utils
import glob
from scipy.misc import imresize

class Series:

    def __init__(self,reference_h5_fn,vidx=None,layer_names=['ISOS','COST']):
        self.tag_template = '%s_%03d'
        self.reference_filename = reference_h5_fn
        self.data_path = os.path.split(reference_h5_fn)[0]
        self.reference_h5 = H5(reference_h5_fn,mode='r')
        self.reference_h5.catalog()
        self.layer_names = layer_names

        self.n_vol = self.reference_h5.get('/config/n_vol').value
        self.n_slow = self.reference_h5.get('/config/n_slow').value
        self.n_fast = self.reference_h5.get('/config/n_fast').value

        print self.n_vol,self.n_slow,self.n_fast
        
        if vidx is None:
            # if the caller doesn't know which volume to use,
            # show all of them and then quit
            stack = np.zeros((self.n_vol,self.n_slow,self.n_fast))
            for vidx in range(self.n_vol):
                for layer_name in layer_names:
                    stack[vidx,:,:] = stack[vidx,:,:]+self.reference_h5['projections'][layer_name][vidx,:,:]/float(len(layer_names))
            clims = np.percentile(stack[:,50:-50,50:-50],(5,99))
            for vidx in range(self.n_vol):
                plt.figure()
                test = stack[vidx,:,:]
                plt.imshow(test,cmap='gray',interpolation='none',aspect='auto',clim=clims)
                plt.colorbar()
                plt.title('volume %d'%vidx)
            plt.show()
            sys.exit()
                        

        self.reference_tag = self.tag_template%(os.path.split(self.reference_filename)[1].replace('.hdf5',''),vidx)
        h5fn = self.tag_template%(self.reference_filename.replace('.hdf5',''),vidx)+'_series.hdf5'
        self.h5 = H5(h5fn)

        self.n_vol = self.reference_h5['config']['n_vol'].value
        self.n_fast = self.reference_h5['config']['n_fast'].value
        self.n_slow = self.reference_h5['config']['n_slow'].value

        if not self.h5.has(self.reference_tag):
            # write 0's and 1's for shifts and goodnesses, respectively, for the reference
            self.h5.put('/%s/x_shifts'%self.reference_tag,np.zeros(self.n_slow))
            self.h5.put('/%s/y_shifts'%self.reference_tag,np.zeros(self.n_slow))
            self.h5.put('/%s/goodnesses'%self.reference_tag,np.ones(self.n_slow))
            self.h5.put('/%s/reference'%self.reference_tag,[1])

        self.reference_vidx = vidx
        
        stack = np.zeros((len(layer_names),self.n_slow,self.n_fast))
        for idx,layer_name in enumerate(self.layer_names):
            stack[idx,:,:] = self.reference_h5['projections'][layer_name][self.reference_vidx,:,:]

        self.reference = np.mean(stack,axis=0)
        self.reference_h5.close()


    def add(self,filename,vidx,slowmin=None,slowmax=None,fastmin=None,fastmax=None,overwrite=False,oversample_factor=5,strip_width=3.0,do_plot=False):
        
        print 'Adding %s, volume %d.'%(filename,vidx)
        if slowmin is None:
            slowmin = 0
        if slowmax is None:
            slowmax = self.n_slow

        if fastmin is None:
            fastmin = 0
        if fastmax is None:
            fastmax = self.n_fast
            
        target_tag = self.tag_template%(os.path.split(filename)[1].replace('.hdf5',''),vidx)

        if self.h5.has(target_tag) and not overwrite:
            print 'Series already has entry for %s.'%target_tag
            return

        target = self.get_image(filename,vidx)
        reference = self.reference
        y,x,g = utils.strip_register(target[slowmin:slowmax,fastmin:fastmax],reference[slowmin:slowmax,fastmin:fastmax],oversample_factor,strip_width,do_plot=do_plot)
        
        self.h5.put('/%s/x_shifts'%target_tag,x)
        self.h5.put('/%s/y_shifts'%target_tag,y)
        self.h5.put('/%s/goodnesses'%target_tag,g)
        self.h5.put('/%s/reference'%target_tag,[0])

    def get_image(self,filename,vidx):
        target_h5 = H5(filename)
        stack = np.zeros((len(self.layer_names),self.n_slow,self.n_fast))
        for idx,layer_name in enumerate(self.layer_names):
            stack[idx,:,:] = target_h5['projections'][layer_name][vidx,:,:]
        return np.mean(stack,axis=0)

    def get_image_tag(self,tag):
        fn,vidx = self.tag_to_filename_vidx(tag)
        return self.get_image(fn,vidx)

    def imshow_tag(self,tag,fig=None):
        im = self.get_image_tag(tag)
        self.imshow(im)


        
    def imshow(self,im,fig=None):
        if fig is None:
            plt.figure()
        else:
            plt.clf()
            plt.cla()
        plt.imshow(im,cmap='gray',clim=np.percentile(im,(5,99.5)))
    
    def catalog(self):
        self.h5.catalog()

    def close(self):
        self.h5.close()


    def fix_yshifts(self):
        sys.exit('Are you sure you want to do this? You probably already fixed it!')
        keys = self.h5.keys()
        for k in keys:
            yshifts = self.h5['%s/y_shifts'%k][:]
            new_yshifts = yshifts - np.arange(self.n_slow)
            self.h5.put('%s/y_shifts_with_offset'%k,yshifts)
            del self.h5.h5['%s/y_shifts'%k]
            self.h5.put('%s/y_shifts'%k,new_yshifts)
        sys.exit()

    def render(self,goodness_threshold=None,correlation_threshold=None,slowmin=None,slowmax=None,fastmin=None,fastmax=None,overwrite=False,keylist=None,oversample_factor=5,do_plot=False):

        if goodness_threshold is None:
            # goodness_threshold = np.median(self.get_all_goodnesses())
            goodness_threshold = 0.0

        if correlation_threshold is None:
            correlation_threshold = 0.0

        if slowmin is None:
            slowmin = 0
        if slowmax is None:
            slowmax = self.n_slow

        if fastmin is None:
            fastmin = 0
        if fastmax is None:
            fastmax = self.n_fast

        if keylist is None:
            keys = self.h5.keys()
        else:
            keys = keylist

        sign = -1
            
        # first, find the minimum and maximum x and y shifts
        xmin = np.inf
        xmax = -np.inf
        ymin = np.inf
        ymax = -np.inf
        for k in keys:
            goodnesses = self.h5['%s/goodnesses'%k][:]
            xshifts = sign*self.h5['%s/x_shifts'%k][:]
            yshifts = sign*self.h5['%s/y_shifts'%k][:]

            xshifts = np.squeeze(xshifts)
            yshifts = np.squeeze(yshifts)

            xshifts,yshifts,goodnesses = self.filter_registration(xshifts,yshifts,goodnesses)

            yshifts = yshifts + np.arange(self.n_slow)
            
            valid = np.where(goodnesses>=goodness_threshold)[0]

            if len(valid):
                xshifts = xshifts[valid]
                yshifts = yshifts[valid]

                newxmin = np.min(xshifts)
                newxmax = np.max(xshifts)
                newymin = np.min(yshifts)
                newymax = np.max(yshifts)
                
                xmin = min(xmin,newxmin)
                xmax = max(xmax,newxmax)
                ymin = min(ymin,newymin)
                ymax = max(ymax,newymax)


        xmin = np.round(xmin*oversample_factor)
        xmax = np.round(xmax*oversample_factor)
        dx = xmax-xmin
        xoffset = xmin
        width = self.n_fast*oversample_factor + dx

        ymin = np.round(ymin*oversample_factor)
        ymax = np.round(ymax*oversample_factor)
        dy = ymax-ymin
        yoffset = ymin
        height = oversample_factor + dy

        sum_image = np.zeros((height,width))
        counter_image = np.zeros((height,width))

        ref_oversampled = imresize(self.reference,int(round(oversample_factor*100)),interp='nearest')

        x1 = round(sign*xoffset)
        x2 = x1+ref_oversampled.shape[1]
        y1 = round(sign*yoffset)
        y2 = y1+ref_oversampled.shape[0]
        fig = plt.figure()
        all_corr_coefs = []
        
        for k in keys:
            print k
            goodnesses = self.h5['%s/goodnesses'%k][:]
            xshifts = sign*self.h5['%s/x_shifts'%k][:]
            yshifts = sign*self.h5['%s/y_shifts'%k][:]
            xshifts = np.squeeze(xshifts)
            yshifts = np.squeeze(yshifts)

            xshifts,yshifts,goodnesses = self.filter_registration(xshifts,yshifts,goodnesses)

            valid = np.where(goodnesses>=goodness_threshold)[0]

            im = self.get_image_tag(k)

            if (not any(xshifts)) and (not any(yshifts)):
                block = imresize(im,int(oversample_factor*100),interp='nearest')
                bsy,bsx = block.shape
                x1 = -xoffset
                y1 = -yoffset
                x2 = x1+bsx
                y2 = y1+bsy
                sum_image[y1:y2,x1:x2] = sum_image[y1:y2,x1:x2] + block
                counter_image[y1:y2,x1:x2] = counter_image[y1:y2,x1:x2] + 1.0

            else:

                if len(valid):
                    for v in valid:
                        line = im[v,fastmin:fastmax]
                        line = np.expand_dims(line,0)
                        block = imresize(line,int(oversample_factor*100),interp='bilinear')
                        bsy,bsx = block.shape
                        x1 = np.round(xshifts[v]*oversample_factor-xoffset)
                        y1 = v*oversample_factor+np.round(yshifts[v]*oversample_factor-yoffset)
                        x2 = x1+bsx
                        y2 = y1+bsy

                        im1 = sum_image[y1:y2,x1:x2].ravel()
                        if np.min(im1)==np.max(im1):
                            corr = 1.0
                        else:
                            corr_valid = np.where(im1)[0]
                            corr = np.corrcoef(im1[corr_valid],block.ravel()[corr_valid])[1,0]
                        all_corr_coefs.append(corr)

                        if corr>correlation_threshold:
                            sum_image[y1:y2,x1:x2] = sum_image[y1:y2,x1:x2] + block
                            counter_image[y1:y2,x1:x2] = counter_image[y1:y2,x1:x2] + 1.0

            if do_plot:
                temp = counter_image.copy()
                temp[np.where(temp==0)] = 1.0
                av = sum_image/temp
                plt.clf()
                plt.axes([0,0,.6,1.0])
                plt.cla()
                plt.imshow(av,cmap='gray',clim=np.percentile(av,(1,99.5)),interpolation='none')
                plt.colorbar()
                plt.axes([.6,.5,.4,.4])
                plt.cla()
                plt.imshow(counter_image)
                plt.colorbar()
                plt.axes([.7,.1,.25,.35])
                plt.cla()
                try:
                    plt.hist(all_corr_coefs,100)
                except:
                    pass
                #plt.colorbar()
                plt.pause(.0000000001)


        plt.figure()
        plt.subplot(1,3,1)
        plt.imshow(ref_oversampled,cmap='gray',clim=np.percentile(ref_oversampled,(1,99.5)),interpolation='none')
        temp = counter_image.copy()
        temp[np.where(temp==0)] = 1.0
        av = sum_image/temp
        plt.subplot(1,3,2)
        plt.imshow(av,cmap='gray',clim=np.percentile(av,(5,99.5)),interpolation='none')
        plt.subplot(1,3,3)
        plt.imshow(counter_image)
        plt.colorbar()

        plt.savefig('%s_rendered.png'%self.reference_tag,dpi=300)
        plt.show()
        

    def filter_registration(self,xshifts,yshifts,goodnesses,xmax=10,ymax=10):

        
        xvalid = np.abs(xshifts - np.median(xshifts))<=xmax
        yvalid = np.abs(yshifts - np.median(yshifts))<=ymax
        invalid = np.logical_or(1-xvalid,1-yvalid)


        # for xv,yv,iv in zip(xvalid,yvalid,invalid):
        #     print xv,yv,iv
        # sys.exit()
        
        # plt.plot(xshifts)
        # plt.plot(yshifts)
        # plt.plot((goodnesses-np.mean(goodnesses))*10000)
        # plt.show()
        # sys.exit()

        goodnesses[invalid] = -np.inf

        return xshifts,yshifts,goodnesses
        
    def goodness_histogram(self):
        all_goodnesses = self.get_all_goodnesses()
        plt.hist(all_goodnesses,100)
        plt.show()
        
    def get_all_goodnesses(self,exclude_reference=True):
        all_goodnesses = []
        keys = self.h5.keys()
        for k in keys:
            goodnesses = self.h5['%s/goodnesses'%k][:]
            if np.prod(goodnesses)==1.0 and exclude_reference:
                continue
            all_goodnesses.append(goodnesses)
        return np.array(all_goodnesses).ravel()


    def tag_to_filename_vidx(self,k):
        last_underscore = k.rfind('_')
        fn = k[:last_underscore]+'.hdf5'
        full_filename = os.path.join(self.data_path,fn)
        vidx = int(k[last_underscore+1:])
        return full_filename,vidx
    
    def show_images(self):

        keys = self.h5.keys()
        for k in keys:
            plt.cla()
            im = self.get_image(full_filename,vidx)
            self.imshow(im)
            plt.title(k)
            plt.pause(.1)
        plt.close()
        
    def find_reference(self):
        keys = self.h5.keys()
        for k in keys:
            if self.h5['%s/reference'%k].value:
                return k

if __name__=='__main__':

    ref_fn = '/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.11.21_cones/14_13_13-1T_500.hdf5' # volume 7

    s = Series(ref_fn,vidx=7,layer_names=['ISOS','COST'])
    #s.show_images()
    #s.goodness_histogram()

    keylist = ['14_13_13-1T_500_%03d'%7]
    for k in range(12):
        if not k==7:
            keylist.append('14_13_13-1T_500_%03d'%k)
    print keylist
    s.render(keylist=None,fastmin=30)
    s.close()
    sys.exit()
    files = glob.glob('/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.11.21_cones/*.hdf5')
    for f in files:
        for vidx in range(12):
            s.add(f,vidx=vidx,fastmin=20)


    s.catalog()
    sys.exit()
    fn = '/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.11.21_cones/14_01_03-1T.hdf5' # volume 0
    tfn = '/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.11.21_cones/14_19_25-1T_500_0.hdf5' # volume 0

    
    fn2 = '/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.08.16/14_44_47-4T_1000.hdf5' # volume 8
    #fn = '/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.08.16/14_44_08-4T_500.hdf5' # volume
    
    s = Series(fn,vidx=1,layer_names=['ISOS'])
    s.add(fn,vidx=0,fastmin=20)
    s.catalog()
    s.close()
    
