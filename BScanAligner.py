from octopod import *
import os,logging,sys
logging.basicConfig(level='INFO')


# THIS MODULE IS FAR FROM FINISHED!!! IT'S A MESS RIGHT NOW.

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
        self.h5.delete(sv_key)
        try:
            zvec = self.h5.get(sv_key)[:]
        except Exception as e:
            self.logger.info('z_align_rough: no %s to delete from h5. This is expected.'%sv_key)
            self.logger.info('z_align_rough: Generating.')
            
            proj = np.mean(np.abs(self.h5.get('/processed_data')[iVol,:,:,:]),axis=1)
            plt.imshow(proj)
            plt.show()
            proj = proj/np.std(proj)
            sy,sx = proj.shape
            refs = []
            shift_vectors = []
            xc_vectors = []
            
            for k in range(0,sy,sy/5):
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

        nv,ns,nd,nf = self.h5.get_shape('/processed_data')
        

        
        
        # slow1 = iVol*self.nSlow
        # slow2 = slow1 + self.nSlow
        # zoff = np.max(zvec)
        
        # outvol = np.zeros((ns,nf,nd+zoff),dtype=np.complex64)
        # for k in range(ns):
        #     outvol[k,:,zvec[k]:zvec[k]+nd] = self.h5['processed'][iVol,k,:,:]


        # try:
        #     zvec = self.h5['/fast_axis_alignment_vector']
        # except Exception as e:
        #     try:
        #         del self.h5['/fast_axis_alignment_vector']
        #     except Exception as e:
        #         pass
        #     self.logger.info('No fast axis alignment vector in h5. Generating.')

        #     proj = np.mean(np.abs(outvol),axis=0)
        #     proj = proj/np.std(proj)
        #     sy,sx = proj.shape
        #     refs = []
        #     shift_vectors = []
        #     xc_vectors = []
            
        #     for k in range(sy/5,sy,sy/5):
        #         refs.append(proj[k,:])
                
        #     for ref in refs:
        #         projf = np.fft.fft(proj,n=sx*2,axis=1)
        #         reff = np.fft.fft(ref,n=sx*2)
        #         xc = np.abs(np.fft.fftshift(np.fft.ifft(projf*np.conj(reff),axis=1),axes=1))
        #         shift_vec = np.argmax(xc,axis=1)
        #         shift_vec = shift_vec - np.median(shift_vec)
                
        #         shift_vectors.append(shift_vec)
        #         xc_vectors.append(np.max(xc,axis=1))

        #     shift_vectors = np.array(shift_vectors)
        #     xc_vectors = np.array(xc_vectors)

        #     shift_vector = np.median(shift_vectors,axis=0)
        #     xc_vector = np.median(xc_vectors,axis=0)

        #     shift_vector = -shift_vector
        #     zvec = shift_vector - np.min(shift_vector)
        #     self.h5['/fast_axis_alignment_vector'] = zvec


        # ns,nf,nd = outvol.shape
        # zoff = np.max(zvec)
        
        # outvol2 = np.zeros((ns,nf,nd+zoff),dtype=np.complex64)
        # for k in range(nf):
        #     outvol2[:,k,zvec[k]:zvec[k]+nd] = outvol[:,k,:]

        # return outvol2



if __name__=='__main__':

    h5 = H5('./oct_test_volume/oct_test_volume_2T.hdf5')
    bsa = BScanAligner(h5)
    bsa.z_align_rough()
