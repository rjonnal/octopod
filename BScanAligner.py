from octopod import *
import os,logging,sys
logging.basicConfig(level='INFO')


class BScanAligner:

    def __init__(self,H5):
        self.h5 = H5
        self.logger = logging.getLogger(__name__)
        self.logger.info('Creating BScanAligner object.')

    def z_align_rough(self,iVol=0,do_plot=True):
        # this method returns a coarsely aligned volume
        # if '/bscan_alignment_vector' exists, it uses this to
        # align the b-scans; if not, it creates that vector first.
        sv_key = '/axial_alignment/slow_vector'
        fv_key = '/axial_alignment/fast_vector'
        ff_key = '/axial_alignment/fast_flattened_volume'
        sf_key = '/axial_alignment/slow_flattened_volume'
        #self.h5.delete(sv_key)
        #self.h5.delete(fv_key)
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


        # for k in range(3):
        #     plt.figure()
        #     plt.imshow(np.mean(np.abs(outvol),axis=k),interpolation='none',aspect='auto')

        # plt.show()


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
            outvol2 = self.h5.get(sf_key)[:]
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

                
            self.h5.put(sf_key,outvol2)

            
        for k in range(3):
            plt.figure()
            plt.imshow(np.mean(np.abs(outvol2),axis=k))
        plt.show()



if __name__=='__main__':

    #h5 = H5('./oct_test_volume/oct_test_volume_2T.hdf5')
    h5 = H5('/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.04.12_2_old/IMG_RE_4TR_10.hdf5')
    bsa = BScanAligner(h5)
    bsa.z_align_rough()
