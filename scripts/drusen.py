from octopod import *
from matplotlib import pyplot as plt
import sys,os
import logging
logging.basicConfig(level=logging.INFO)
from octopod.utils import nxcorr2

fm = FileManager()

start_index = 50 # which image to use as base, among 100 images in each dataset

def test_nxcorr2():
    target = np.random.rand(100,100)
    template = target[20:30,20:30]
    xshift,yshift,peakVal,nxcval = nxcorr2(template,target)
    print xshift,yshift,peakVal,nxcval


def label_druse_across_sets(subject_name):
    files = fm.get('hroct')
    files = [f for f in files if f.find(subject_name)>-1]

    druse_id = None
    for fn in files:
        d = Dataset(fn)
        h5 = d.get_h5_handle()
        im = h5['processed_data'][0][start_index]
        fig = plt.figure()

        h5.require_group('drusen')
        points = []
        def on_click(event):
            print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(event.button, event.x, event.y, event.xdata, event.ydata)
            points.append([event.xdata,event.ydata])

        fig.canvas.mpl_connect('button_press_event',on_click)

        ih = plt.imshow(np.log(np.abs(im)))
        ih.set_clim((9,12))
        plt.colorbar()
        plt.title(os.path.split(fn)[-1])
        plt.show()
        if druse_id is None:
            druse_id = raw_input('Please enter an ID for this druse: ')
        try:
            del h5['drusen'][druse_id]
        except Exception as e:
            pass
        h5['drusen'].create_group(druse_id)
        h5['drusen'][druse_id].create_dataset('coordinates',data=points)

def move_coordinates_in_h5(subject_name):
    files = fm.get('hroct')
    files = [f for f in files if f.find(subject_name)>-1]
    druse_id = None
    for fn in files:
        d = Dataset(fn)
        h5 = d.get_h5_handle()
        keys = h5['drusen'].keys()
        for key in keys:
            arr = np.zeros((3,2))
            arr[:] = h5['drusen'][key][:]
            del h5['drusen'][key]
            h5['drusen'].create_group(key)
            h5['drusen'][key].create_dataset('coordinates',data=arr)

def translation(im0, im1):
    """Return translation vector to register images."""
    shape = im0.shape
    f0 = np.fft.fft2(im0)
    f1 = np.fft.fft2(im1)
    ir = abs(np.fft.ifft2((f0 * f1.conjugate()) / (abs(f0) * abs(f1))))
    goodness = np.max(ir)
    t0, t1 = np.unravel_index(np.argmax(ir), shape)
    if t0 > shape[0] // 2:
        t0 -= shape[0]
    if t1 > shape[1] // 2:
        t1 -= shape[1]
    return [t0, t1, goodness]

def track_druse_through_set(subject_name,druse_id,border=20):
    files = fm.get('hroct')
    files = [f for f in files if f.find(subject_name)>-1]

    for fn in files:
        d = Dataset(fn)
        h5 = d.get_h5_handle()
        print h5['drusen'].keys()
        coords_ref = h5['drusen'][druse_id]['coordinates']
        xs = coords_ref[:,0]
        ys = coords_ref[:,1]
        x1 = round(np.min(xs)-border)
        x2 = round(np.max(xs)+border)
        y1 = round(np.min(ys)-border)
        y2 = round(np.max(ys)+border)
        temp = h5['processed_data'][0][start_index]
        template = np.zeros(temp.shape)
        #template[y1:y2,x1:x2] = np.abs(temp[y1:y2,x1:x2])
        template[:,:] = np.abs(temp[:,:])
        for k in range(50,100):
            target = np.abs(h5['processed_data'][0][k])
            print translation(template,target)

def register_images(subject_name):
    pass
register_images('Carmen','2p5deg_temporal')
#track_druse_through_set('Carmen','2p5deg_temporal')


    
#move_coordinates_in_h5('Carmen')
#label_druse_across_sets('Carmen')
