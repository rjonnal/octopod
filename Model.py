import sys,os
import logging
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from utils import translation,translation1,autotrim_bscan,find_peaks,shear
from octopod.Misc import H5
import logging
import octopod_config as ocfg
logging.basicConfig(level=logging.DEBUG)


if False:
    x = np.random.rand(100,100)
    a = x[5:,:]
    b = x[:-5,:]
    print translation(a,b)

class Model:
    """A model of gross retinal reflectance, used for segmenting
    and labeling A-scans."""

    def __init__(self,h5,debug=False):
        """Initialize model. May pass an h5py.File object or .hdf5 filename.
        The file or object must contain a 'processed_data' dataset containing
        at least one volume."""
        self.logger = logging.getLogger(__name__)
        if type(h5)==str:
            self.h5 = H5(h5)
            self.logger.info('Opening file %s.'%self.h5)
        else:
            self.h5 = h5
        try:
            self.profile = self.h5.get('model/profile')[:]
        except Exception as e:
            self.logger.info('Model does not exist in h5 file.')
            self.profile = self.make_model(debug=debug)

    def blur(self,bscan,kernel_width=5):
        return sp.signal.convolve2d(bscan,np.ones((1,kernel_width)),mode='same')/float(kernel_width)
        
    def make_model(self,vidx=0,debug=False):
        self.logger.info('Making model...')
        avol = np.abs(self.h5.get('processed_data')[vidx,:,:,:])
        nSlow,nFast,nDepth = avol.shape
        
        if False:
            # make a test volume with very clean translations.
            avol = np.ones(avol.shape)
            zpos = 100
            rad = 3
            zpos_vec = []
            for iSlow in range(nSlow):
                avol[iSlow,zpos-rad:zpos+rad,:] = avol[iSlow,zpos-rad:zpos+rad,:] + 100.0
                zpos = np.abs(np.round(zpos + np.random.randn()))
                zpos_vec.append(zpos)

                
        # crop to avoid edge effects in the registration:
        if avol.shape[2]>10:
            avol = avol[:,:,5:-5]

        tx_vec = [0.0]
        ty_vec = [0.0]
        g_vec = [1.0]

        template = None

        # compute an intensity threshold for the first template image,
        # to avoid using 

        thresh = np.mean(avol)*.75

        for iSlow in range(1,nSlow):
            last = self.blur(avol[iSlow-1,:,:])

            bright_enough = np.mean(last)>=thresh
            
            if debug:
                self.logger.debug('Frame %d of %d. Sufficiently bright: %d'%(iSlow,nSlow-1,bright_enough))
            
            if template is None:
                if bright_enough:
                    template = last
                    counter = np.ones(template.shape)
                    sy,sx = template.shape
                else:
                    continue
            
            target = self.blur(avol[iSlow,:,:])
            
            if False:
                plt.subplot(1,2,1)
                plt.cla()
                plt.imshow(template,interpolation='none',aspect='auto')
                #plt.colorbar()
                plt.subplot(1,2,2)
                plt.cla()
                plt.imshow(target,interpolation='none',aspect='auto')
                #plt.colorbar()
                plt.pause(.1)
            
            tx,ty,g = translation(target,template,debug=False)
            if debug:
                self.logger.debug('x:%d; y:%d; goodness:%0.3f'%(tx,ty,g))
                
            # ty is the amount to shift target to align with template
            # let's say template's profile is [0 1 0 0] and target's is
            # [0 0 1 0]. ty = -1 in this case.
            if ty<0:
                put0 = -ty
                put1 = sy
                get0 = 0
                get1 = ty
            elif ty>0:
                put0 = 0
                put1 = -ty
                get0 = ty
                get1 = sy
            else:
                get0 = 0
                get1 = sy
                put0 = 0
                put1 = sy

            template[put0:put1,:] = template[put0:put1,:] + target[get0:get1,:]
            counter[put0:put1,:] = counter[put0:put1] + 1.0

            if debug:
                plt.cla()
                plt.plot(np.mean(template/counter,axis=1))
                plt.pause(.1)


        template = template/counter
        template = shear(template,2)

        model_profile = np.mean(template,axis=1)
        
        self.write_profile(model_profile)
        return model_profile

    def write_profile(self,model_profile):
        self.h5.require_group('model')
        self.h5.put('model/profile',model_profile)

    def click_crop(self):
        global clicks
        clicks = []
        
        fig = plt.figure(figsize=(10,6))
        
        def onclick(event):
            global clicks
            newclick = round(event.xdata)
            clicks.append(newclick)
            plt.axvline(newclick)
            plt.draw()
            
        cid = fig.canvas.mpl_connect('button_press_event',onclick)
        
        plt.plot(self.profile)
        plt.show()

        if len(clicks)>=2:

            z1 = np.min(clicks)
            z2 = np.max(clicks)
            self.profile[:z1] = 0.0
            self.profile[z2+1:] = 0.0
            self.write_profile(self.profile)
            


    def get_label_dict(self):
        out = {}
        try:
            for key in self.h5.get('model/labels').keys():
                out[key] = self.h5.get('model/labels')[key]
        except Exception as e:
            print e
        return out

    def clear_labels(self):
        self.h5.delete('model/labels')


    def get_volume_labels(self,z_tolerance=2):
        nvol,nslow,ndepth,nfast = self.h5.get('processed_data').shape
        offset_matrix = self.h5.get('model/z_offsets')
        goodness_matrix = self.h5.get('model/z_offset_goodness')

        label_keys = self.h5.get('model/labels').keys()
        labels = {}
        volume_labels = {}
        for key in label_keys:
            labels[key] = self.h5.get('model/labels')[key].value
            volume_labels[key] = np.zeros((nvol,nslow,nfast))
            
        for ivol in range(nvol):
            avol = np.abs(self.h5.get('processed_data')[ivol,:,:,:])
            for islow in range(nslow):
                for ifast in range(nfast):
                    test = avol[islow,:,ifast]
                    offset = offset_matrix[ivol,islow,ifast]
                    goodness = goodness_matrix[ivol,islow,ifast]
                    for key in label_keys:
                        model_z_index = labels[key]
                        temp = np.zeros(test.shape)
                        temp[...] = test[...]
                        temp[:(model_z_index-offset-z_tolerance)] = 0.0
                        temp[(model_z_index-offset+z_tolerance)+1:] = 0.0
                        volume_labels[key][ivol,islow,ifast] = np.argmax(temp)

        for key in label_keys:
            location = 'model/volume_labels/%s'%key
            self.h5.put(location,volume_labels[key])
            plt.figure()
            plt.imshow(volume_labels[key][0,:,:],interpolation='none',aspect='auto')
            plt.title(key)
            plt.colorbar()

        plt.figure()
        plt.imshow(offset_matrix[0,:,:],interpolation='none',aspect='auto')
        plt.title('offset')
        plt.colorbar()
        plt.figure()
        plt.imshow(goodness_matrix[0,:,:],interpolation='none',aspect='auto')
        plt.title('goodness')
        plt.colorbar()
        plt.show()
        
                
    def write_axial_alignment(self):
        nvol,nslow,ndepth,nfast = self.h5.get('processed_data').shape
        offset_matrix = self.h5.make('model/z_offsets',(nvol,nslow,nfast),dtype='i2')
        goodness_matrix = self.h5.make('model/z_offset_goodness',(nvol,nslow,nfast),dtype='f8')
        om,gm = self.align_volumes()
        offset_matrix[...] = om[...]
        goodness_matrix[...] = gm[...]

    def align_volumes(self):
        nvol,nslow,ndepth,nfast = self.h5.get('processed_data').shape
        offset_matrix = np.zeros((nvol,nslow,nfast))
        goodness_matrix = np.zeros((nvol,nslow,nfast))
        
        for ivol in range(nvol):
            offset,goodness = self.align_volume(vidx=ivol)
            offset_matrix[ivol,:,:] = offset
            goodness_matrix[ivol,:,:] = goodness

        return offset_matrix,goodness_matrix
        
    def align_volume(self,vidx=0,rad=5):
        avol = np.abs(self.h5.get('processed_data')[vidx,:,:,:])
        avol = np.swapaxes(avol,0,1)
        ndepth,nslow,nfast = avol.shape
        FF,SS = np.meshgrid(np.arange(nfast),np.arange(nslow))
        offset_submatrix = np.zeros((nslow,nfast))
        goodness_submatrix = np.zeros((nslow,nfast))
        profile = self.profile
        if len(profile)>ndepth:
            profile = profile[:ndepth]
        if ndepth>len(profile):
            avol = avol[:len(profile),:,:]
            ndepth,nslow,nfast = avol.shape
        
        for islow in range(nslow):
            self.logger.info('B-scan %d in volume %d.'%(islow,vidx))
            for ifast in range(nfast):
                self.logger.info('A-scan %d.'%ifast)
                ff = FF - ifast
                ss = SS - islow
                d = np.sqrt(ff**2 + ss**2)
                d[np.where(d>rad)]=0
                d[np.where(d)]=1
                d = d/np.sum(d)
                fvol = avol*d
                test = np.mean(np.mean(fvol,axis=2),axis=1)
                offset,goodness = translation1(profile,test,debug=False)
                offset_submatrix[islow,ifast] = offset
                goodness_submatrix[islow,ifast] = goodness
        return offset_submatrix,goodness_submatrix
                

    def click_label(self,smoothing=1):
        if smoothing>1:
            working_profile = sp.signal.convolve(self.profile,np.ones((smoothing)),mode='same')/float(smoothing)
        else:
            working_profile = self.profile
        
        # find peaks and troughs:
        gthresh = 1.0/smoothing

        peaks = list(find_peaks(working_profile,gradient_threshold=gthresh))+list(find_peaks(-working_profile,gradient_threshold=gthresh))
        peaks = sorted(peaks)
        
        idx = 0
        z = np.arange(len(working_profile))
        done = False or not len(peaks)


        fig = plt.figure(figsize=(20,6))
        for key in plt.rcParams.keys():
            if key[:6]=='keymap':
                plt.rcParams[key] = ''
        
        global current_x,current_label,label_dict
        label_dict = self.get_label_dict()
        current_x = 0
        current_label = ''

        def plot_at(x):
            global current_label

            plt.subplot(131)
            plt.cla()
            plt.imshow(np.tile(working_profile,(512,1)).T,interpolation='none',aspect='auto')
            
            plt.subplot(132)
            plt.cla()
            plt.plot(z,working_profile)
            plt.plot(z[x],working_profile[x],'ks')
            valid = np.where(working_profile)[0]
            plt.xlim((valid[0],valid[-1]))
            plt.autoscale(False)
            for label in label_dict.keys():
                label_z = z[label_dict[label]]
                plt.text(label_z,working_profile[label_z],label,ha='center',va='bottom')

            plt.subplot(133)
            plt.cla()
            plt.plot(z,working_profile)
            plt.plot(z[x],working_profile[x],'ks')
            plt.xlim((z[x]-10,z[x]+10))
            z1 = max(0,x-10)
            z2 = min(len(z),x+10)
            ymin = np.min(working_profile[z1:z2])
            ymax = np.max(working_profile[z1:z2])
            plt.ylim((ymin,ymax))
            plt.title(current_label)
            plt.draw()


        def onclick(event):
            global current_x
            # print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
            #     event.button, event.x, event.y, event.xdata, event.ydata)
            current_x = round(event.xdata)
            plot_at(current_x)

        def onpress(event):
            global current_x,current_label,label_dict
            if event.key=='right':
                current_x = (current_x + 1)%len(working_profile)
            elif event.key=='ctrl+right':
                current_x = (current_x + 20)%len(working_profile)
            elif event.key=='left':
                current_x = (current_x - 1)%len(working_profile)
            elif event.key=='ctrl+left':
                current_x = (current_x - 20)%len(working_profile)
            elif event.key=='shift':
                pass
            elif event.key=='/':
                pass
            elif event.key in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ':
                current_label = current_label + event.key.upper()
            elif event.key=='backspace':
                current_label = current_label[:-1]
            elif current_label=='delete':
                label_dict = {}
            elif event.key=='enter':
                label_dict[current_label] = current_x
                print label_dict
                current_label = ''
            plot_at(current_x)
            

        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        pid = fig.canvas.mpl_connect('key_press_event', onpress)
        
        plot_at(current_x)
        
        plt.show()

        self.h5.require_group('model')
        self.h5.require_group('model/labels')
        for key in label_dict.keys():
            self.h5.put('model/labels/%s'%key,label_dict[key])
            
        mdb = H5(ocfg.model_database)
        did = self.h5.get('IDs/dataset_id').value
        # did is the primary key for the model, but we'll also save eccentricity
        did_key = '%d'%did
        mdb.require_group(did_key)

        si = self.h5.get('eccentricity/superior_inferior').value
        nt = self.h5.get('eccentricity/nasal_temporal').value
        
        radial_distance = np.sqrt(si**2+nt**2)

        mdb.put('%s/superior_inferior'%did_key,si)
        mdb.put('%s/nasal_temporal'%did_key,nt)
        mdb.put('%s/radial_distance'%did_key,radial_distance)
        mdb.put('%s/profile'%did_key,self.profile)

        mdb.require_group('%s/labels'%did_key)
        for key in label_dict.keys():
            mdb.put('%s/labels/%s'%(did_key,key),label_dict[key])
        


def test():
    h5 = H5('./oct_test_volume/oct_test_volume_2T.hdf5')
    m = Model(h5,True)

    #m.clear_labels()
    #m.click_crop()
    #m.click_label()
    #m.write_axial_alignment()
    m.get_volume_labels()
    
if __name__=='__main__':
    test()
