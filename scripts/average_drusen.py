from octopod import *
import os


fm = FileManager()
files = fm.get('hroct',['Carmen'])


def save_to_png(fn,image,clims,scale=1.0,boxx=[],boxy=[],dpi=150):
    sy,sx = image.shape
    sy_inches = float(sy)/float(dpi)*scale
    sx_inches = float(sx)/float(dpi)*scale
    fig = plt.figure(figsize=(sx_inches,sy_inches))
    ax = plt.axes([0,0,1,1])
    imh = ax.imshow(image)
    imh.set_clim(clims)
    if len(boxx):
        plt.autoscale(False)
        plt.plot(boxx,boxy,'y-',lw=1)
    fig.savefig(fn)
    

for f in files:

    wd_name = os.path.splitext(f)[0]+'_working'
    
    d = Dataset(f)
    h5 = d.get_h5_handle()

    drusen_ids = h5['drusen'].keys()

    for druse_id in drusen_ids:

        avifn = os.path.join(wd_name,'%s.avi'%druse_id)
        avgfn = os.path.join(wd_name,'%s_average.png'%druse_id)
        bscanfn = os.path.join(wd_name,'%s_bscan.png'%druse_id)

        track = h5['drusen'][druse_id]['tracking'][:]
        ntrack,ncol = track.shape
        idx,x1,x2,y1,y2,goodness = track[0,:]
        im_sum = h5['processed_data'][0][idx].copy()
        im_sum = im_sum[y1:y2,x1:x2]
        im_sum = np.abs(im_sum)

        border = 5
        boxx = [x1-border,x2+border,x2+border,x1-border,x1-border]
        
        boxy = [y1-border,y1-border,y2+border,y2+border,y1-border]

        reference_index = h5['drusen'][druse_id]['reference_index'][()]
        bscan = np.log(np.abs(h5['processed_data'][0,reference_index,:,:]))

        save_to_png(bscanfn,bscan,(9,12),scale=1.0,boxx=boxx,boxy=boxy)
        print bscanfn
        continue
        
        for t in range(1,ntrack):
            print '%s: %d of %d'%(os.path.split(f)[1],t+1,ntrack)
            idx,x1,x2,y1,y2,goodness = track[t,:]
            im = h5['processed_data'][0][idx].copy()
            im = im[y1:y2,x1:x2]
            im = np.abs(im)
            if im_sum.shape[1]==im.shape[1] and im_sum.shape[0]==im.shape[0]:
                im_sum = im_sum + im
            plt.imshow(im,interpolation='none')
            plt.pause(.1)
            plt.cla()
        plt.close()
        plt.imshow(im_sum,interpolation='none')
        plt.show()

