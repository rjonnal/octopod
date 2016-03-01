from octopod import *
from matplotlib import pyplot as plt
import sys,os
import logging
logging.basicConfig(level=logging.INFO)
from octopod.utils import translation

fm = FileManager()

start_index = 40 # which image to use as base, among 100 images in each dataset

def test_nxcorr2():
    target = np.random.rand(100,100)
    template = target[20:30,20:30]
    xshift,yshift,peakVal,nxcval = nxcorr2(template,target)
    print xshift,yshift,peakVal,nxcval


def label_druse_across_sets(subject_name):
    files = fm.get('hroct')
    files = [f for f in files if f.find(subject_name)>-1]
    
    druse_id = None
    for fn in files[::-1]:
        d = Dataset(fn)
        h5 = d.get_h5_handle()

        index_to_use = start_index - 1
        chosen = False

        while not chosen:
            index_to_use = index_to_use + 1
            im = h5['processed_data'][0][index_to_use]
            ih = plt.imshow(np.log(np.abs(im)))
            ih.set_clim((9,12))
            tp = np.arange(0,1000,100)
            tv = ['%0.1f'%x for x in tp/1000.0*5.0/0.3-3.0]
            plt.xticks(tp)
            plt.gca().set_xticklabels(tv)
            plt.pause(.1)
            chosen = raw_input('Use [%s]? '%druse_id)=='y'

        
        plt.close()
        im = h5['processed_data'][0][index_to_use]
        fig = plt.figure()

        h5.require_group('drusen')
        points = []
        def on_click(event):
            print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(event.button, event.x, event.y, event.xdata, event.ydata)
            points.append([event.xdata,event.ydata])

        fig.canvas.mpl_connect('button_press_event',on_click)

        ih = plt.imshow(np.log(np.abs(im)))
        ih.set_clim((9,12))
        plt.xticks(tp)
        plt.gca().set_xticklabels(tv)
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
        h5['drusen'][druse_id].create_dataset('reference_index',data=[index_to_use])
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

def track_druse_through_all_sets(subject_name,druse_id,border=20):
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
        template[:,x1:x2] = np.abs(temp[:,x1:x2])
        #template[:,:] = np.abs(temp[:,:])
        out = []
        for k in range(100):
            target = np.abs(h5['processed_data'][0][k])
            xshift,yshift,goodness = translation(template,target)
            d = np.sqrt(xshift**2+yshift**2)
            print k,d,goodness
            myx1 = x1 - xshift
            myx2 = x2 - xshift
            myy1 = y1 - yshift
            myy2 = y2 - yshift
            debug = False
            if debug:
                plt.subplot(2,2,1)
                plt.cla()
                plt.imshow(template[y1:y2,x1:x2])
                plt.subplot(2,2,2)
                plt.cla()
                plt.imshow(target[myy1:myy2,myx1:myx2])
                plt.subplot(2,2,3)
                plt.cla()
                plt.imshow(template)
                plt.subplot(2,2,4)
                plt.cla()
                plt.imshow(target)
                plt.pause(.1)
            if goodness>.02:
                out.append([k,myx1,myx2,myy1,myy2,goodness])
                print '\t added'
        try:
            del h5['drusen'][druse_id]['tracking']
        except:
            pass
        h5['drusen'][druse_id].create_dataset('tracking',data=np.array(out))

def track_druse_through_set(filename,druse_id,border=20):
    d = Dataset(filename)
    h5 = d.get_h5_handle()
    print h5['drusen'].keys()
    coords_ref = h5['drusen'][druse_id]['coordinates']
    xs = coords_ref[:,0]
    ys = coords_ref[:,1]
    x1 = round(np.min(xs)-border)
    x2 = round(np.max(xs)+border)
    y1 = round(np.min(ys)-border)
    y2 = round(np.max(ys)+border)
    refidx = h5['drusen'][druse_id]['reference_index'][0]
    print refidx
    temp = h5['processed_data'][0][refidx]
    template = np.zeros(temp.shape)
    template[:,x1:x2] = np.abs(temp[:,x1:x2])
    out = []
    for k in range(100):
        target = np.abs(h5['processed_data'][0][k])
        xshift,yshift,goodness = translation(template,target)
        d = np.sqrt(xshift**2+yshift**2)
        print k,d,goodness
        myx1 = x1 - xshift
        myx2 = x2 - xshift
        myy1 = y1 - yshift
        myy2 = y2 - yshift
        debug = False
        if debug:
            plt.subplot(2,2,1)
            plt.cla()
            plt.imshow(template[y1:y2,x1:x2])
            plt.subplot(2,2,2)
            plt.cla()
            plt.imshow(target[myy1:myy2,myx1:myx2])
            plt.subplot(2,2,3)
            plt.cla()
            plt.imshow(template)
            plt.subplot(2,2,4)
            plt.cla()
            plt.imshow(target)
            plt.pause(.1)
        out.append([k,myx1,myx2,myy1,myy2,goodness])
    try:
        del h5['drusen'][druse_id]['tracking']
    except:
        pass
    h5['drusen'][druse_id].create_dataset('tracking',data=np.array(out))
        
        
label_druse_across_sets('Edric')

carmen_ids = ['1_deg_nasal','0p5_deg_temporal','2_deg_temporal','5_deg_temporal']
edric_ids = ['0p5_deg_nasal','5_deg_temporal']

files = fm.get('hroct',['Carmen'])
for f in files:
    for did in carmen_ids:
        track_druse_through_set(f,did)

