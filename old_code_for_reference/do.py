from OCTTools import *
from time import sleep
from multiprocessing import Pool
from time import sleep
from matplotlib import pyplot as plt
import sys,os

skip_server = False
search_string = 'McDermott_Kyle'

s = Searcher()
for x in s.unprocessed+s.processed:
    if x.find(search_string)>-1:
        print x
        
if True:# basic post-processing off campus (no server access)
    for f in s.unprocessed:
        if f.find(search_string)>-1:
            p = Processor(f,skip_server)
            p.optimizeDispersionCoefficients(plot=False)
            p.archiveProcessed(plot=False)
            p.makeMovie()
            
if True:# b-scan alignment
    for f in s.processed:
        if f.find(search_string)>-1:
            p = Processor(f,skip_server)
            temp = os.path.split(os.path.split(p.diagnosticPngDirectory)[0])
            temp0 = temp[0]
            temp1 = temp[1]
            first = os.path.split(temp0)[1]
            second = temp1
            outfn = '%s_%s.png'%(first,second)
            # try:
            #     del p.h5['/bscan_alignment_vector']
            # except Exception as e:
            #     print e
            # continue
            vol = p.z_align_rough()
            vol = np.abs(vol)
            ns,nf,nd = vol.shape
            prof = np.mean(np.mean(vol,axis=1),axis=0)

            left = prof[:-2]
            center = prof[1:-1]
            right = prof[2:]
            peaks = np.where(np.logical_and(center>right,center>left))[0]+1
            peak_vals = prof[peaks]
            valid = np.where(peak_vals>2000)[0]
            pmax = peaks[valid[0]]

            
            func = np.mean
            vol = vol[:,:,pmax-3:pmax+4]

            plt.figure(figsize=(8,3))
            plt.subplot(1,2,1)
            plt.plot(prof)
            plt.plot(peaks,prof[peaks],'k.')
            plt.axvspan(pmax-3,pmax+3)
            plt.subplot(1,2,2)
            plt.imshow(func(vol,axis=2))
            #plt.axis('equal')
            plt.axis('tight')
            plt.colorbar()
            outfn2 = os.path.join(p.diagnosticPngDirectory,'en_face.png')
            plt.savefig(outfn,dpi=300)
            plt.savefig(outfn2,dpi=300)
            plt.pause(.1)
            #plt.close()
            continue
            for k in zrange:
                plt.cla()
                plt.imshow(np.abs(vol[:,:,k]))
                plt.pause(.0001)

if False:# output PNG stack
    for f in s.processed:
        if f.find('Werner_Eli')>-1:
            p = Processor(f,True)

#p = Processor('H:/AOOCT_Data/2014/2014.08.26/Francis_Dennis/LE/1mmLineNFL_300NR_000IR.unp')
#p.process()
#p.makeMovie()
#sys.exit()
#p.optimizeDispersionCoefficients(plot=False)
#p.writeGlobalDispersionCoefs((1,2,3))
#sys.exit()
#s.cleanup()
#s.makeConfigFiles(doAll=True)
#s.manualView(doAll=True)
#s.archivePairwiseShifts()


sys.exit()

p.process()
p.makeMovie()

#p.deleteDiagnostics(False)
#p.deleteLogs(False)
#p.optimizeDispersionCoefficients()
#p.showDispersionCompensation()






sys.exit()


#p = Processor('H:/Share/afocal_aooct_data/Data/2013.12.20/Jonnal_Ravi/RE/5.0N/20um-ON---0.unp')

#data = p.getFrame(0,0)
#data = p.process(data,cropDC=True)
#sys.exit()


def pfunk(pid):
    sleep(.5)
    p = Processor('/dose/Share/afocal_aooct_data/Data/2013.12.20/Jonnal_Ravi/RE/5.0N/20um-ON---0.unp')
    #p.deleteDiagnostics()
    p.dispersionCoefs = [0.0,0.0,0.0,0.0]
    return p.optimizeDispersionCoefficientsDebug(plot=False)

pool = Pool(processes=6)
result = pool.map(pfunk,range(6))
output = []
for item in result:
    dc3 = item[0][0]
    dc2 = item[0][1]
    preSNR = item[1]
    postSNR = item[2]
    output.append([dc3,dc2,preSNR,postSNR])

output = np.array(output)
np.savetxt('results.txt',output)

plt.figure()
plt.plot(output[:,0],output[:,3]/output[:,2],'ko')
plt.figure()
plt.plot(output[:,1],output[:,3]/output[:,2],'bo')
plt.show()

sys.exit()
