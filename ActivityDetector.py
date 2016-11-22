from octopod import *
import os,logging,sys
logging.basicConfig(level='INFO')


class ActivityDetector:

    def __init__(self,h5):
        self.h5 = h5
        self.logger = logging.getLogger(__name__)
        self.logger.info('Creating ActivityDetector object.')
        self.n_vol = self.h5['config']['n_vol'].value
        self.n_slow = self.h5['config']['n_slow'].value
        self.n_fast = self.h5['config']['n_fast'].value
        self.n_depth = self.h5['config']['n_depth'].value

        
    def project_series(self):
        fast_project = []
        depth_project = []

        slow_project_started = False
        for i_vol in range(self.n_vol):
            print i_vol
            for i_slow in range(self.n_slow):
                b = self.h5['processed_data'][i_vol,i_slow,:,:]
                ab = np.abs(b)
                if not slow_project_started:
                    slow_project = ab
                    slow_project_started = True
                    slow_count = 0
                else:
                    slow_project = slow_project + ab
                    slow_count = slow_count + 1
                

                fast_project.append(np.mean(ab,axis=1))
                depth_project.append(np.mean(ab,axis=0))
                
        fast_project = np.array(fast_project)
        depth_project = np.array(depth_project)
        slow_project = slow_project/float(slow_count)

        plt.figure()
        plt.imshow(depth_project,interpolation='none',aspect='auto')
        plt.colorbar()
        plt.figure()
        plt.imshow(fast_project,interpolation='none',aspect='auto')
        plt.colorbar()
        plt.figure()
        plt.imshow(slow_project,interpolation='none',aspect='auto')
        plt.colorbar()
        plt.show()
        
        
    def align_volumes(self,do_plot=False,keep_intermediate=False):
        nv,ns,nd,nf = self.h5.get('processed_data').shape
        for iVol in range(nv):
            self.align_volume(iVol,do_plot=do_plot)

        depths = []
        for iVol in range(nv):
            shape = self.h5['/axial_alignment/flattened_volume_%03d'%iVol].shape
            depths.append(shape[1])

        self.logger.info('Depths of flattened volumes: %s.'%(','.join(['%d'%d for d in depths])))
        
        new_nd = np.max(depths)
        self.logger.info('Max depth: %d'%new_nd)
        
        flat_block = np.zeros((nv,ns,new_nd,nf),dtype=np.complex64)
        self.logger.info('Creating %dx%dx%dx%d volume.'%(nv,ns,new_nd,nf))
        
        for iVol in range(nv):
            flat_vol = self.h5['/axial_alignment/flattened_volume_%03d'%iVol][:]
            self.logger.info('Shape of volume to insert: %s'%(','.join(['%d'%d for d in flat_vol.shape])))
            dmax = flat_vol.shape[1]
            self.logger.info('Inserting at %d,:,:%d,:.'%(iVol,dmax))
            flat_block[iVol,:,:dmax,:] = flat_vol
        self.h5.put('/flattened_data',flat_block)
        if not keep_intermediate:
            self.h5.delete('axial_alignment')
            self.h5.repack()
        
            
    def align_volume(self,iVol=0,do_plot=False):
        # this method returns a coarsely aligned volume
        # if '/bscan_alignment_vector' exists, it uses this to
        # align the b-scans; if not, it creates that vector first.
        sv_key = '/axial_alignment/slow_vector_%03d'%iVol
        fv_key = '/axial_alignment/fast_vector_%03d'%iVol
        ff_key = '/axial_alignment/fast_flattened_volume_%03d'%iVol
        f_key = '/axial_alignment/flattened_volume_%03d'%iVol
        #self.h5.delete(sv_key)
        #self.h5.delete(fv_key)
        #self.h5.delete(ff_key)
        #self.h5.delete(f_key)

        try:
            zvec = self.h5.get(fv_key)[:]
        except Exception as e:
            self.logger.info('z_align_rough: no %s to delete from h5. This is expected.'%fv_key)
            self.logger.info('z_align_rough: Generating.')

            # first we average in the slow dimension
            proj = np.mean(np.abs(self.h5.get('/processed_data')[iVol,:,:,:]),axis=0).T
            proj = proj/np.std(proj)
            sy,sx = proj.shape
            refs = []
            shift_vectors = []
            xc_vectors = []


            chunk_size = sy/5
            for k in range(0,sy,sy/5):
                print k+chunk_size/2
                refs.append(proj[k+chunk_size/2,:])

            for ref in refs:
                projf = np.fft.fft(proj,n=sx*2,axis=1)
                reff = np.fft.fft(ref,n=sx*2)
                xc = np.abs(np.fft.fftshift(np.fft.ifft(projf*np.conj(reff),axis=1),axes=1))
                shift_vec = np.argmax(xc,axis=1)
                shift_vec = shift_vec - np.median(shift_vec)
                
                shift_vectors.append(shift_vec)
                xc_vectors.append(np.max(xc,axis=1))

            shift_vectors = np.array(shift_vectors)
            xc_vectors = np.array(xc_vectors)

            
            shift_vector = np.median(shift_vectors,axis=0)
            xc_vector = np.median(xc_vectors,axis=0)

            shift_vector = -shift_vector
            zvec = shift_vector - np.min(shift_vector)
            self.h5.put(fv_key,zvec)

        try:
            outvol = self.h5.get(ff_key)[:]
        except Exception as e:
            self.logger.info('%s'%e)
            nv,ns,nd,nf = self.h5.get_shape('/processed_data')
        
            zoff = np.max(zvec)

            # build a new volume where each fast scan is flattened
            self.logger.info('z_align_rough: Fetching volume from h5.')
            invol = self.h5.get('/processed_data')[iVol,:,:,:]
            outvol = np.zeros((ns,nd+zoff,nf),dtype=np.complex64)
            self.logger.info('z_align_rough: Aligning slow-scans to one another.')
            for k in range(nf):
                if k%20==0:
                    self.logger.info('z_align_rough: %d percent done.'%(float(k)/float(nf)*100))
                outvol[:,zvec[k]:zvec[k]+nd,k] = invol[:,:,k]
            self.h5.put(ff_key,outvol)


        try:
            zvec = self.h5.get(sv_key)[:]
        except Exception as e:
            self.logger.info('No %s in h5. Generating.'%sv_key)
            
            proj = np.mean(np.abs(outvol),axis=2)
            proj = proj/np.std(proj)
            
            sy,sx = proj.shape
            refs = []
            shift_vectors = []
            xc_vectors = []
            
            for k in range(sy/5,sy,sy/5):
                refs.append(proj[k,:])
                
            for ref in refs:
                projf = np.fft.fft(proj,n=sx*2,axis=1)
                reff = np.fft.fft(ref,n=sx*2)
                xc = np.abs(np.fft.fftshift(np.fft.ifft(projf*np.conj(reff),axis=1),axes=1))
                shift_vec = np.argmax(xc,axis=1)
                shift_vec = shift_vec - np.median(shift_vec)
                
                shift_vectors.append(shift_vec)
                xc_vectors.append(np.max(xc,axis=1))

            shift_vectors = np.array(shift_vectors)
            xc_vectors = np.array(xc_vectors)

            shift_vector = np.median(shift_vectors,axis=0)
            xc_vector = np.median(xc_vectors,axis=0)

            shift_vector = -shift_vector
            zvec = shift_vector - np.min(shift_vector)
            self.h5.put(sv_key,zvec)
            
        try:
            outvol2 = self.h5.get(f_key)[:]
        except Exception as e:
            self.logger.info('%s'%e)
            self.logger.info('z_align_rough: Assembling volume with flattened slow scans.')
            ns,nd,nf = outvol.shape
            zoff = np.max(zvec)

            outvol2 = np.zeros((ns,nd+zoff,nf),dtype=np.complex64)
            for k in range(ns):
                if k%20==0:
                    self.logger.info('z_align_rough: %d percent done'%(float(k)/float(ns)*100))
                outvol2[k,zvec[k]:zvec[k]+nd,:] = outvol[k,:,:]
            self.h5.put(f_key,outvol2)


        if do_plot:
            for k in range(3):
                plt.figure()
                plt.imshow(np.mean(np.abs(outvol2),axis=k))
            plt.show()



if __name__=='__main__':

    #h5 = H5('./oct_test_volume/oct_test_volume_2T.hdf5')
    h5 = H5('/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.04.12_2_old/IMG_RE_4TR_10.hdf5')
    bsa = BScanAligner(h5)
    bsa.z_align_rough()
