from tools import *
from scipy.cluster.vq import whiten,kmeans2
from scipy import signal
import sys
from multiprocessing import Pool
from time import time
import os
from matplotlib import pyplot as plt
from time import sleep
import png
from scipy.misc import imresize

class GrowableMatrix:

    def __init__(self):
        self._data = np.empty((1,1))
        self._x1 = 0
        self._x2 = 1
        self._y1 = 0
        self._y2 = 1


    def addRow(self,row,y,x):
        #print 'data shape:',self._data.shape
        if x<self._x1:
            sy,sx = self._data.shape
            # grow matrix to the left
            dx = self._x1 - x
            temp = np.zeros((sy,sx+dx))
            self._x1 = self._x1 - dx
            temp[:,dx:dx+self._data.shape[1]] = self._data
            self._data = temp

        if x+len(row)>self._x2:
            sy,sx = self._data.shape
            # grow matrix to the right
            dx = (x+len(row)) - self._x2
            temp = np.zeros((sy,sx+dx))
            #print 'case 2 temp shape:',temp.shape
            self._x2 = self._x2 + dx
            temp[:,:self._data.shape[1]] = self._data
            self._data = temp

        if y<self._y1:
            sy,sx = self._data.shape
            dy = self._y1 - y
            temp = np.zeros((sy+dy,sx))
            self._y1 = self._y1 - dy
            temp[dy:dy+self._data.shape[0],:] = self._data
            self._data = temp

        if y>=self._y2:
            sy,sx = self._data.shape
            dy = y - self._y2 + 1
            temp = np.zeros((sy+dy,sx))

            #print 'case 4 temp shape:',temp.shape
            self._y2 = self._y2 + dy
            #print 'case 4 data shape:',self._data.shape

            temp[:self._data.shape[0],:] = self._data
            self._data = temp

        self.genericInsert(row,y,x)


    def genericInsert(self,row,y,x):
        yinsert = y - self._y1
        xinsert = x - self._x1
        xinsert2 = xinsert + len(row)
        self._data[yinsert,xinsert:xinsert+len(row)] = self._data[yinsert,xinsert:xinsert+len(row)] + row

    def getData(self):
        return self._data


    def replace(self,oldval,newval):
        self._data[np.where(self._data==oldval)]=newval

class Averager:
    def __init__(self,label):
        self._sum = GrowableMatrix()
        self._count = GrowableMatrix()
        self._label = label
        self._layers = []

    def add(self,row,y,x,z=0):
        counter = np.ones(row.shape)
        self._sum.addRow(row,y,x)
        self._count.addRow(counter,y,x)
        if z<len(self._layers):
            self._layers[z].addRow(row,y,x)
        elif z==len(self._layers):
            self._layers.append(GrowableMatrix())
            self._layers[z].addRow(row,y,x)
        else:
            sys.exit('Averager.add: z value (%d) exceeds length of self._layers list (%d).'%(z,len(self._layers)))
            

    def getLayer(self,idx):
        return self._layers[idx]

    def get(self,frac=0.0):

        self._count.replace(0,1)

        s = self._sum._data
        c = self._count._data

        cnew = np.zeros(c.shape)
        cnew[:] = c[:]

        cnew[np.where(cnew==0)]=1
            
        im = s.astype(np.float64)/cnew.astype(np.float64)
#        im = self.gaussianBlur(im,15,0.7)
        return im,c


    def gaussianBlur(self,im,N,sigma):
        nRange = np.linspace(-(N-1)/2,(N-1)/2,N)
        xx,yy = np.meshgrid(nRange,nRange)
        d2 = xx**2+yy**2
        g = np.exp(-(d2/(2*sigma**2)))
        g = g/np.sum(g)
        return signal.convolve2d(im,g,'same')


    def getLabel(self):
        return self._label

    def show(self,frac=0.0):
        # plt.subplot(131)
        # implt = plt.imshow(self._sum._data,aspect='auto',interpolation='none')
        # implt.set_cmap('gray')
        # plt.colorbar()

        # plt.subplot(132)
        # implt = plt.imshow(self._count._data,aspect='auto',interpolation='none')
        # implt.set_cmap('gray')
        # plt.colorbar()
        
        #plt.subplot(133)
        im = self.get(frac)
        immax = im[np.where(im)].max()*0.9
        immin = im[np.where(im)].min()*1.1
        implt = plt.imshow(self.get(),interpolation='none',aspect='auto')
        implt.set_cmap('gray')
        implt.set_clim((immin,immax))
        plt.colorbar()
        #plt.show()
        

class VolumeRegistrar:

    def __init__(self,referenceSet,volumeIndex,regLayerIndex,regLayerOffsets,layersToAverage,layersToAverageLabels,nxcThreshold=0.5):
        

        self._layersToAverage = layersToAverage
        self._layersToAverageLabels = layersToAverageLabels

        refFnList = glob.glob(os.path.join(referenceSet._processingCache,'layerIntensity_%d_%d*.npy'%(regLayerIndex,volumeIndex)))
        refFn = refFnList[0]

        refStack = np.load(refFn)
        self._ref = refStack[regLayerOffsets,:,:].mean(0)


        self._outputDirectory = os.path.join(referenceSet.getWorkingDirectory(),'vrOutput')

        if not os.path.exists(self._outputDirectory):
            os.makedirs(self._outputDirectory)

        self._regLayerOffsets = regLayerOffsets
        self._regLayerIndex = regLayerIndex
        self._nxcThreshold = nxcThreshold
        self._tags = []
        self._setListFilename = os.path.join(self._outputDirectory,'setList.txt')
        self._fileListFilename = os.path.join(self._outputDirectory,'fileList.txt')


        self._hasData = os.path.exists(self._fileListFilename)

        self._averagers = []
        for lta,ltal in zip(self._layersToAverage,self._layersToAverageLabels):
            self._averagers.append(Averager(ltal))

        self._dataSetCount = 0

        self._singles = None # this will later be replaced by a stack of single images

    def add(self,ds,smooth=False):

        self._dataSetCount = self._dataSetCount + 1
        if not self._hasData:
            fid = open(self._fileListFilename,'w')
            fid.write('')
            fid.close()

            fid = open(self._setListFilename,'w')
            fid.write('')
            fid.close()
            


        self._tags.append(ds.getTag())

        fid = open(self._fileListFilename,'a')
        fid.write(ds.getFilename()+'\n')
        fid.close()


        layerDirectory = ds._processingCache
        

        nVol = ds.getNVol()
        nLayersToAverage = len(self._layersToAverage)

        for iVol in range(nVol):
            for iLayer in range(nLayersToAverage):
                tarFnList = glob.glob(os.path.join(layerDirectory,'layerIntensity_%d_%d*.npy'%(iLayer,iVol)))
                tarFn = tarFnList[0]
                temp = np.load(tarFn)[self._regLayerOffsets,:,:].mean(0) # load the images, and average over the bright layers

                if smooth:
                    temp = ndimage.gaussian_filter(temp,.75)
                
                if iLayer==0:
                    sy,sx = temp.shape
                    tarStack=np.zeros((nLayersToAverage,sy,sx))
                
                tarStack[iLayer,:,:] = temp

            if self._singles is None:
                self._singles = tarStack


            self._addVolume(tarStack,ds.getTag(),iVol)


    def _addVolume(self,newImageStack,tag,volumeIndex,forceRebuild=False,interactive=False):

        target = newImageStack[self._regLayerIndex,:,:]
        sy,sx = target.shape

        peakFnEnd = 'peaks_%s_%04d.npy'%(tag,volumeIndex)
        shiftFnEnd = 'shifts_%s_%04d.npy'%(tag,volumeIndex)
        peakFn = os.path.join(self._outputDirectory,peakFnEnd)
        shiftFn = os.path.join(self._outputDirectory,shiftFnEnd)


        cacheExists = os.path.exists(peakFn) and os.path.exists(shiftFn)
        if cacheExists:
            print 'Cached peaks and shifts exist in %s and %s. Loading from these.'%(peakFn,shiftFn)


        try:
            peaksTemp = np.load(peakFn)
            shiftsTemp = np.load(shiftFn)
        except Exception as e:
            print e
            peaksTemp = np.zeros((sy,sy))
            shiftsTemp = np.zeros((sy,sy))

        plt.ion()
        plt.clf()
        plt.subplot(131)
        plt.imshow(peaksTemp)
        plt.colorbar()

        plt.subplot(132)
        plt.imshow(target)
        plt.colorbar()

        plt.subplot(133)
        plt.imshow(self._ref)
        plt.colorbar()
        plt.draw()

        
        lIdxTag = '%04d'%self._regLayerIndex


        fid = open(self._setListFilename,'a')
        fid.write('%s_%04d_%s\n'%(tag,volumeIndex,lIdxTag))
        fid.close()


        y = True
        n = False

        if interactive:
            doThisPair = input('Compute this pair? (y/n) ')
        else:
            doThisPair = True

        if not doThisPair:
            shifts = np.ones((sy,sy))*np.inf
            peaks = np.ones((sy,sy))*np.inf*(-1)
            #np.save(peakFn,peaks)
            #np.save(shiftFn,shifts)
            return

        plt.ioff()

        if (not cacheExists) or forceRebuild:

            shifts = np.zeros((sy,sy),dtype=np.float16)
            peaks = np.zeros((sy,sy),dtype=np.float16)

            for tidx in range(sy):
                print 'Line %d of %d'%(tidx+1,sy)
                tarLine = target[tidx,:]
                for ridx in range(sy):
#                    print 'Line %d of %d'%(tidx+1,sy),' - %d of %d'%(ridx+1,sy)
                    refLine = self._ref[ridx,:]
                    try:
                        shift,peakVal = self.nxcorr1(tarLine,refLine,True)
                    except Exception as e:
                        print e,'setting shift and peak to nan'
                        shift = np.nan
                        peakVal = np.nan
                    #print shift,peakVal

                    shifts[tidx,ridx] = np.float16(shift)
                    peaks[tidx,ridx] = np.float16(peakVal)

                if not (tidx+1)%10:
                    plt.ion()
                    plt.subplot(121)
                    plt.imshow(peaks)
                    plt.subplot(122)
                    plt.imshow(shifts.astype(np.float16)-shiftsTemp)
                    plt.draw()

            np.save(peakFn,peaks)
            np.save(shiftFn,shifts)
            self._hasData = True
        else:
            peaks = np.load(peakFn).astype(np.float16)
            shifts = np.load(shiftFn).astype(np.float16)
            
        info = self.conditionShifts(peaks,shifts)


        infoarr = np.array(info)
        np.savetxt('./volume_%02d_map.txt'%volumeIndex,infoarr)


        for idx,layer in enumerate(self._layersToAverage):
            for tidx,tup in enumerate(info):
                if tidx==0:
                    lastX = 0
                    lastY = 0
                tarLine = newImageStack[layer,tup[0],:]
                x = tup[2]
                y = tup[1]
                d = np.sqrt((x-lastX)**2+(y-lastY)**2)
                if d<=5:
                    self._averagers[idx].add(tarLine,y,x,volumeIndex)
                lastX = x
                lastY = y


    def bmpScale(self,im,count,scaleFactor):

        if count is None:
            count = np.ones(im.shape)

        immax = im[np.where(count)].max()
        immin = im[np.where(count)].min()

        im = (im - immin)/(immax - immin)
        im = im.clip(0,1)
        im = np.round(im*255)
        im = np.uint8(im)

        if scaleFactor - 1:
            im = imresize(im,scaleFactor,'nearest','L')

        return im



    def makeImages(self,scaleFactor=1.5,border=5):

        nLayers = len(self._layersToAverage)

        singles = self._singles
        
        for idx,layer in enumerate(self._layersToAverage):
            
            label = self._averagers[idx].getLabel()
            im,count = self._averagers[idx].get()
            

            outfn = os.path.join(self._outputDirectory,'%s_avg.npy'%label)
            np.save(outfn,im)


            
            if idx==1 or idx==2:
                tempoutdir = os.path.join('/','home','rjonnal','Dropbox','figures','src','averaged_layers')
                singleImage = singles[idx,:,:]
                averageImage = im

                np.save('singleImage_%d.npy'%idx,singleImage)
                np.save('averageImage_%d.npy'%idx,averageImage)
                

                
                if False:
                    cropSize = 100

                    axcrop1 = 196
                    axcrop2 = axcrop1+100
                    aycrop1 = 153
                    aycrop2 = aycrop1+100

                    sxcrop1 = 24
                    sxcrop2 = sxcrop1+100
                    sycrop1 = 147
                    sycrop2 = sycrop1+100

                    subAverage = averageImage[aycrop1:aycrop2,axcrop1:axcrop2]
                    subSingle = singleImage[sycrop1:sycrop2,sxcrop1:sxcrop2]

                    anxc = np.abs(np.fft.fftshift(np.fft.fft2(subAverage)))
                    anxcl = np.log(anxc)
                    plt.figure()
                    imh = plt.imshow(anxcl)
                    # plt.colorbar()
                    imh.set_clim((6,14))

                    snxc = np.abs(np.fft.fftshift(np.fft.fft2(subSingle)))
                    snxcl = np.log(snxc)
                    plt.figure()
                    imh = plt.imshow(snxcl)
                    # plt.colorbar()
                    imh.set_clim((6,14))

                    plt.figure()
                    plt.imshow(np.log(anxc/snxc))
                    plt.colorbar()


                    plt.show()


                    sys.exit()
            

            im = self.bmpScale(im,count,scaleFactor)
            sim = self.bmpScale(singles[idx,:,:],None,scaleFactor)

            if idx==0:
                aHeight,aWidth = im.shape
                sHeight,sWidth = sim.shape
                height = max(aHeight,sHeight)
                canvas = np.zeros((height*nLayers+border*(nLayers-1),aWidth+sWidth+border),dtype=np.uint8)

            canvas[idx*(border+aHeight):idx*(border+aHeight)+aHeight,:aWidth] = im
            canvas[idx*(border+aHeight):idx*(border+aHeight)+sHeight,aWidth+border:] = sim

            if idx==2:
                for lidx,layer in enumerate(self._averagers[idx]._layers):
                    
                    layer_image = layer.getData()
                    layer_image[np.where(layer_image<1)] = 0
                    lsy,lsx = layer_image.shape
                    ar = float(lsy)/float(lsx)
                    plt.clf()
                    plt.axes([0,0,1,1])
                    imh = plt.imshow(layer_image)
                    plt.xticks([])
                    plt.yticks([])
                    plt.axis('equal')
                    #plt.axis('tight')
                    immin = layer_image[np.where(layer_image)].min()
                    layer_image[np.where(layer_image==0)] = immin
                    llim = immin*1.5
                    ulim = layer_image.max()*.67
                    imh.set_clim((llim,ulim))
                    plt.savefig('/home/rjonnal/Dropbox/meetings/2014/ARVO/poster/figures/cone_mosaic_corrected_%d.png'%lidx,dpi=600,facecolor='k',edgecolor='none')

                    plt.show()



        maxCount = count.max()
        fn = os.path.join(self._outputDirectory,'averaged_layers_%d.png'%maxCount)
        png.from_array(canvas,'L').save(fn)
        
        cfn = os.path.join(self._outputDirectory,'averaged_layers_count_%d.png'%maxCount)
        plt.clf()
        plt.imshow(count)
        plt.colorbar()
        plt.savefig(cfn)
        

    def doWeirdThings(self):
        setList = []
        fid = open(self._setListFilename)
        setList = fid.readlines()
        fid.close()

        peakList = []
        shiftList = []
        for set in setList:
            items = set.split('_')
            tag = ''
            for substr in items[:-3]:
                tag = tag + substr + '_'
            tag = tag + items[-3] + '_' + items[-2] + '.npy'
            peakFn = os.path.join(self._outputDirectory,'peaks_'+tag)
            shiftFn = os.path.join(self._outputDirectory,'shifts_'+tag)
            peakList.append(peakFn)
            shiftList.append(shiftFn)

        for peakFn,shiftFn in zip(peakList,shiftList):
            print peakFn,shiftFn
            peaks = np.load(peakFn)
            shifts = np.load(shiftFn)
            self.conditionShifts(peaks,shifts)

    def conditionShifts(self,peaks,shifts,doPlots=False):
        # 5, 2
        # 5, 1/2 good
        # 5, 0/3/4 bad
        kSize = 5
        minD = 1
        kernel = (-1)*np.ones((kSize,kSize))
        for j in range(kSize):
            for k in range(kSize):
                if np.abs(j-k)<=minD:
                    kernel[j,k] = 1.0

        cPeaks = signal.convolve2d(peaks,kernel,'same')/np.sum(kernel)
        cPeaks = signal.convolve2d(cPeaks,np.ones((kSize,kSize)),'same')/(kSize**2)

        if doPlots:
            plt.subplot(221)
            plt.imshow(cPeaks)
            plt.colorbar()
            plt.subplot(222)
            plt.imshow(peaks)
            plt.colorbar()


        pMaxIdxVec = np.argmax(peaks,axis=0)
        cpMaxIdxVec = np.argmax(cPeaks,axis=0)

        pMaxVec = np.max(peaks,axis=0)
        cpMaxVec = np.max(cPeaks,axis=0)

        sVec = np.zeros(pMaxVec.shape)
        for idx,maxIdx in enumerate(cpMaxIdxVec):
            sVec[idx] = shifts[maxIdx,idx]

        valid = np.where(cpMaxVec>self._nxcThreshold)[0]


        if doPlots:
            plt.subplot(223)
            plt.plot(valid,cpMaxIdxVec[valid],'ks')
            plt.subplot(224)
            plt.plot(valid,sVec[valid],'ks')
            plt.show()


        targetReferenceShiftTuples = []
        for idx in valid:
            tarIdx = cpMaxIdxVec[idx]
            refIdx = idx
            shift = shifts[tarIdx,refIdx]
            tup = (tarIdx,refIdx,shift)
            targetReferenceShiftTuples.append(tup)

        return targetReferenceShiftTuples


    def nxcorr1(self,vec1,vec2,doPlots=False):
        '''Returns shift,xc:
        shift is the number of pixels that vec2 must be
          shifted in order to align with vec1
        xc is the cross correlation peak'''

        l1 = len(vec1)
        l2 = len(vec2)

        vec1 = (vec1 - np.mean(vec1))/np.std(vec1)
        vec2 = (vec2 - np.mean(vec2))/np.std(vec2)

        temp1 = np.zeros([l1+l2-1])
        temp2 = np.zeros([l1+l2-1])

        temp1[:l1] = vec1
        temp2[:l2] = vec2

        nxcval = np.real(fftshift(ifft(fft(temp1)*np.conj(fft(temp2)))))

        peakVal = np.max(nxcval)
        peakIdx = np.where(nxcval==peakVal)[0][0]

        if len(nxcval)%2:
            shift = (len(nxcval)-1)/2.0 - peakIdx
        else:
            shift = len(nxcval)/2.0 - peakIdx


        peakVal = peakVal/len(vec1)
        
        if peakVal>1110.6:
            print shift
            plt.subplot(311)
            plt.plot(vec1)
            plt.subplot(312)
            plt.plot(vec2)
            plt.subplot(313)
            plt.plot(nxcval/len(vec1))
            plt.show()

        return shift,peakVal



if __name__=='__main__':
    
    try:
        case = int(sys.argv[1])
    except Exception as e:
        print e,
        print 'Using 0 for case.'
        case = 0

    if case==0:
        datadir = '/dose/Share/afocal_aooct_data/Data/2013.12.20/Jonnal_Ravi/RE/5.0N/'
        
        files = glob.glob(os.path.join(datadir,'20um-ON---0.unp'))

        referenceTag = '20um-ON---0'

        referenceLayerIndex = 2
        referenceLayerOffsets = [8,9,10,11,12]

        referenceVolumeIndex = 4

        layersToAverage=range(4)
        layerNames = ['ELM','ISOS','OST','RPE']


    refset = DataSet(os.path.join(datadir,'%s.unp'%referenceTag))
    
    refsubvol1 = refset.loadLayerSubvolume(1,referenceVolumeIndex)
    refsubvol2 = refset.loadLayerSubvolume(2,referenceVolumeIndex)
    
    subvols1 = []
    subvols2 = []
    mosaics = []
    cxs = []
    cys = []

    for j in range(5):
        cones = refset.loadLayerSubvolume(2,j)[9:12,:,:].mean(0)
        subvols1.append(refset.loadLayerSubvolume(1,j))
        subvols2.append(cones)
        mosaic = Mosaic(cones)
        cy,cx = mosaic.getCoords(doPlots=False)
        cys.append(cy)
        cxs.append(cx)

    refcones = refsubvol2[9:12,:,:].mean(0)
    refMos = Mosaic(refcones)
    

    vr = VolumeRegistrar(refset,referenceVolumeIndex,referenceLayerIndex,referenceLayerOffsets,layersToAverage,layerNames,0.4)

    buildVr = True
    if buildVr:
        for f in files:
            ds = DataSet(f)
            vr.add(ds,smooth=True)
        vr.makeImages()


