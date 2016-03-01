from octopod import *
import os


fm = FileManager()
files = fm.get('hroct',['Carmen'])


def save_to_png(fn,image,clims,scale=1.0,boxx=[],boxy=[],dpi=150):
    fig = layout_for_png(image,clims,scale=scale,boxx=boxx,boxy=boxy,dpi=dpi)
    fig.savefig(fn)
    plt.close(fig)
    
def layout_for_png(image,clims,scale=1.0,boxx=[],boxy=[],dpi=150):
    sy,sx = image.shape
    sy_inches = float(sy)/float(dpi)*scale
    sx_inches = float(sx)/float(dpi)*scale
    fig = plt.figure(figsize=(sx_inches,sy_inches))
    ax = plt.axes([0,0,1,1])
    imh = ax.imshow(image)
    imh.set_clim(clims)
    imh.set_cmap('gray')
    if len(boxx):
        plt.autoscale(False)
        plt.plot(boxx,boxy,'y-',lw=1)
    return fig


clims = [8,13]
for f in files:

    wd_name = os.path.splitext(f)[0]+'_working'
    
    d = Dataset(f)
    h5 = d.get_h5_handle()

    drusen_ids = h5['drusen'].keys()

    for druse_id in drusen_ids:

        avifn = os.path.join(wd_name,'%s.avi'%druse_id)
        avgfn = os.path.join(wd_name,'%s_average.png'%druse_id)
#        if os.path.exists(avgfn):
#            continue
        bscanfn = os.path.join(wd_name,'%s_bscan.png'%druse_id)

        if os.path.exists(avifn):
            print bscanfn+' exists'
            continue
        
        reference_index = h5['drusen'][druse_id]['reference_index'][0]
        track = h5['drusen'][druse_id]['tracking'][:]
        ntrack,ncol = track.shape
        
        idx,x1,x2,y1,y2,goodness = track[reference_index,:]
        
        im_sum = h5['processed_data'][0][idx].copy()
        im_sum = im_sum[y1:y2,x1:x2]
        im_sum = np.log(np.abs(im_sum))

        border = 5
        boxx = [x1-border,x2+border,x2+border,x1-border,x1-border]
        boxy = [y1-border,y1-border,y2+border,y2+border,y1-border]

        bscan = np.log(np.abs(h5['processed_data'][0,reference_index,:,:]))

        save_to_png(bscanfn,bscan,clims,scale=1.0,boxx=boxx,boxy=boxy)
                
        mov = Movie(avifn,cmin=clims[0],cmax=clims[1],make_wmv=False)
        counter = 1.0
        for t in range(1,ntrack):
            print '%s: %d of %d'%(os.path.split(f)[1],t+1,ntrack)
            idx,x1,x2,y1,y2,goodness = track[t,:]
            if goodness>0.017:
                counter = counter + 1
                im = h5['processed_data'][0][idx].copy()
                im = im[y1:y2,x1:x2]
                im = np.log(np.abs(im))
                if im_sum.shape[1]==im.shape[1] and im_sum.shape[0]==im.shape[0]:
                    im_sum = im_sum + im
                fig = layout_for_png(im,clims,scale=5.0)
                try:
                    mov.add(fig)
                except:
                    pass
                plt.close(fig)
                # plt.imshow(im,interpolation='none')
                # plt.pause(.1)
                # plt.cla()
        mov.make()
        im_sum = im_sum/counter
        save_to_png(avgfn,im_sum,clims,scale=5.0)
        #plt.imshow(im_sum,interpolation='none')
        #plt.show()

