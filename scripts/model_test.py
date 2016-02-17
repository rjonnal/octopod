from octopod import Model
import h5py

fn = '/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.02.02/3T_vid_1.hdf5'
#fn = '/home/rjonnal/data/Dropbox/Share/2g_aooct_data/Data/2016.02.02/2T_1.hdf5'

model = Model(fn)
model.make_model(debug=False)
