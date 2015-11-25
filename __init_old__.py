import numpy as np
from matplotlib import pyplot as plt
import sys,os
import time as timemod
import imp
import glob
import h5py
from scipy import interpolate,optimize
from scipy.signal import hilbert,convolve2d
from time import sleep,time
import misc
from movie import Movie
import copy
import fnmatch
import hashlib
import paramiko
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import shutil

try:
    import octtools_config as ocfg
except Exception as e:
    print e
    print 'Attempting to import default configuration.'
    try:
        import octtools_config_default as ocfg
    except Exception as e:
        sys.exit('Exiting: %s.'%e)

class Server:
    def __init__(self,logger=None,dummy=False):
        if logger is None:
            self.log = Logger('./SFTP_log.txt')
        else:
            self.log = logger

        username = ocfg.SERVER_USERNAME
        password = ocfg.SERVER_PASSWORD
        server = ocfg.SERVER_IP
        self.wwwroot = ocfg.SERVER_WWW_ROOT
        port=22
        self.active = False

        if not dummy:
            try:
                self.active = True
                self.log.log('Opening SFTP connection to %s:%d.'%(server,port))
                self.transport = paramiko.Transport((server, port)) 
                self.log.log('Connecting to %s with credentials in octtools_cfg.'%server)
                self.transport.connect(username=username,password=password)
                self.sftp = paramiko.SFTPClient.from_transport(self.transport)

                self.log.log('Opening SSH connection to %s.'%server)
                self.ssh = paramiko.SSHClient()
                self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
                self.ssh.connect(server, username=username,password=password)
            except Exception as e:
                print e
                self.active = False

    def __del__(self):
        print 'Closing SERVER.'
        try:
            self.sftp.close()
        except Exception as e:
            print e

        try:
            self.transport.close()
        except Exception as e:
            print e

    def put(self,src,dest,relative=True):
        retval = 0
        if self.active:

            if relative:
                fulldest = os.path.join(self.wwwroot,dest)
            else:
                fulldest = dest

            self.log.log('Putting %s in %s.'%(src,fulldest))
            try:
                self.sftp.put(src,fulldest)
                retval = 1
            except Exception as e:
                print('FTP.put failed: %s'%e)
        else:
            self.log.log('Server inactive.')

        return retval


    def mkdir(self,dest):
        if self.active:
            fulldest = os.path.join(self.wwwroot,dest)
            self.log.log('Making remote directory %s.'%fulldest)
            try:
                self.sftp.mkdir(fulldest,mode=755)
            except Exception as e:
                print('FTP.mkdir failed: %s'%e)
            self.exec_command('chmod -R 755 %s'%fulldest)


    def exec_command(self,command):
        if self.active:
            self.log.log('Executing command: %s.'%command)
            stdin, stdout, stderr = self.ssh.exec_command(command)
            self.log.log('Server stdout:')
            for line in stdout.readlines():
                self.log.log('\t%s'%line[:-1])

    def addVolume(self,vid,sid,eid,date,eye,info,processed=0):
        if self.active:
            cmd = 'python /home/vsri/www/db/addVolume.py %s %s %s %s %s %s %d'%(vid,sid,eid,date,eye,info,processed)
            self.exec_command(cmd)
            self.mkdir('%s'%sid)
            self.mkdir('%s/%s'%(sid,eid))
            self.mkdir('%s/%s/%s'%(sid,eid,vid))

    def updateProcessed(self,idstr,status):
        if self.active:
            cmd = 'python /home/vsri/www/db/updateProcessed.py %s %d'%(idstr,status)
            self.exec_command(cmd)

    def publish(self):
        if self.active:
            cmd = 'chmod -R 755 /var/www'
            self.exec_command(cmd)

        
            
class ConfigurationTranslator:

    def __init__(self,log=None):
        if log is None:
            self.log = Logger('./log.txt')
        else:
            self.log = log

        self._lookup = {}
        self._lookup['image height'] = 'nFast'
        self._lookup['number of volumes'] = 'nVol'
        self._lookup['number of frames'] = 'nSlow'
        self._lookup['image width'] = 'nDepth'
        self._lookup['xscan amplitude'] = 'xmV'
        self._lookup['yscan volts'] = 'ymV'
        self._lookup['bscanframerate'] = 'bscanRate'
        self._lookup['scan offset'] = 'xOffsetmV'
        self._lookup['yscan offset'] = 'yOffsetmV'

        
    def translate(self,infn,outfn=None):
        if outfn is None:
            outfn = infn.replace('.unp.txt','_cfg.py')

        self.log.log('Translating %s to %s'%(infn,outfn))

        datestamp = timemod.strftime('%Y_%m_%d')
        timestamp = timemod.strftime('%H:%M:%S')
        self._header = '# OCT configuration file automatically generated by \n' \
            '# \tOCTTools.ConfigurationTranslator.translate() \n' \
            '# Time: %s \n' \
            '# Date: %s \n'%(timestamp,datestamp) 

        infid = open(infn,'r')

        outfid = open(outfn,'w')

        outfid.write(self._header)

        line = infid.readline()
        while len(line):
            param,val = line.split(':')
            try:
                outparam = self._lookup[param.lower()]
                val = val.strip()
                outfid.write('%s = %s\n'%(outparam,val))
                self.log.log('\tSetting %s to %s.'%(outparam,val))
            except Exception as e:
                #print '\t',e,
                self.log.log('\t%s could not be parsed'%(line[:-2]))
            line = infid.readline()
            # if we're using this translate method, we must be working with 2G data:
        outfid.write('system_id = \'2G-AO-OCT\'\n')

    def getCfg(self,fn):

        self.path,self.filename = os.path.split(fn)
        self.cfgfn = self.filename.replace('.unp','_cfg.py')
        self.cfgmodule = self.filename.replace('.unp','_cfg')

        # Use Yifan's machine generated config file, translated by
        # ConfigurationTranslater.translate() above, to import
        # a configuration file if none exists.
        yifan_cfgfn = fn + '.txt'

        # try to import a configuration file:
        try:
            self.log.log('Trying to import %s.'%(self.cfgmodule))
            fp,pathname,description = imp.find_module(self.cfgmodule,[self.path])
            cfg = imp.load_module(self.cfgmodule,fp,pathname,description)
            self.log.log('\tsuccess!')
        except Exception as e:
            #print e, '. Attempting to load generic cfg.py.'
            try:
                self.log.log('Trying to import generic cfg.py from %s.'%(self.path))
                fp,pathname,description = imp.find_module('cfg',[self.path])
                cfg = imp.load_module('cfg',fp,pathname,description)
                self.log.log('\tsuccess!')
            except Exception as e:
                #print e, 'Attempting to translate %s.'%yifan_cfgfn
                try:
                    self.log.log('Trying to translate %s into a cfg module.'%yifan_cfgfn)
                    self.translate(yifan_cfgfn)
                    fp,pathname,description = imp.find_module(self.cfgmodule,[self.path])
                    cfg = imp.load_module(self.cfgmodule,fp,pathname,description)
                    self.log.log('\tsuccess!')
                except Exception as e:
                    print e
                    sys.exit('Cannot load any configuration file.')
        return cfg



class Searcher:
    
    def __init__(self):
        droots = ocfg.DATA_ROOTS
        matches = []
        for droot in droots:
            for root, dirnames, filenames in os.walk(droot):
                for filename in fnmatch.filter(filenames, '*.unp'):
                    matches.append(os.path.join(root, filename))

        

        self.unprocessed = []
        self.processed = []
        for fn in matches:
            processed = False

            h5fn = fn.replace('.unp','.hdf5')
            if os.path.exists(h5fn):
                h5 = h5py.File(h5fn)
                try:
                    test = h5['processed']
                    processed = True
                except Exception as e:
                    pass
            if processed:
                self.processed.append(fn)
            else:
                self.unprocessed.append(fn)

        print('Searcher file lists:')
        print('Processed:')
        for fn in self.processed:
            print('\t%s'%fn)
        print('Unprocessed:')
        for fn in self.unprocessed:
            print('\t%s'%fn)

    def cleanup(self):
        todo = self.unprocessed + self.processed
        for fn in todo:
            wdir = fn.replace('.unp','')
            hdf5 = fn.replace('.unp','.hdf5')
            pyc = fn.replace('.unp','_cfg.pyc')
            try:
                print 'removing %s'%wdir
                shutil.rmtree(wdir)
            except Exception as e:
                print e
            
            try:
                print 'removing %s'%hdf5
                os.remove(hdf5)
            except Exception as e:
                print e
                
            try:
                print 'removing %s'%pyc
                os.remove(pyc)
            except Exception as e:
                print e
                
    def optimize(self,doAll=False):
        todo = self.unprocessed
        if doAll:
            todo = todo + self.processed

        for fn in todo:
            try:
                p = Processor(fn)
                p.deleteDiagnostics(False)
                p.deleteLogs(False)
                p.dispersionCoefs = ocfg.DISPERSION_COEFS
                p.optimizeDispersionCoefficients(plot=False)
            except Exception as e:
                print e

    def archiveProcessed(self,doAll=False):
        todo = self.unprocessed
        if doAll:
            todo = todo + self.processed
        
        for fn in todo:
            try:
                plt.ion()
                plt.close('all')
                p = Processor(fn)
                p.archiveProcessed()
            except Exception as e:
                print e
                continue
#            try:
#                p.makeMovie()
#            except Exception as e:
#                print e

    def archivePairwiseShifts(self):
        todo = self.unprocessed + self.processed
        for fn in todo:
            if fn.lower().find('line')>-1:
                try:
                    p = Processor(fn)
                    p.archivePairwiseShifts()
                except Exception as e:
                    print e
                    continue

    def makeMovies(self,doAll=False):
        todo = self.unprocessed + self.processed

        for fn in todo:
            try:
                plt.ion()
                plt.close('all')
                p = Processor(fn)
                p.makeMovie()
            except Exception as e:
                print e
                continue
#            try:
#                p.makeMovie()
#            except Exception as e:
#                print e


    def catalog(self,doAll=False):
        todo = self.unprocessed
        if doAll:
            todo = todo + self.processed
        
        for fn in todo:
            p = Processor(fn)

    def makeConfigFiles(self,doAll=False):
        todo = self.unprocessed
        if doAll:
            todo = todo + self.processed
        
        for fn in todo:
            nVol = 1
            isnfl = fn.upper().find('NFL')>-1
            is200um = fn.upper().find('200UM')>-1
            nlines = os.stat(fn).st_size/2048/2
            nDepth = 2048
            if isnfl and is200um and (nlines/200.0)%1.0==0:
                nSlow = nlines/200
                nFast = 200
            elif isnfl and (nlines/300.0)%1.0==0:
                nSlow = nlines/300
                nFast = 300
            elif (nlines/1000.0)%1.0==0:
                nSlow = nlines/1000
                nFast = 1000
            else:
                nSlow = 1
                nFast = nlines
            print fn,nVol,nSlow,nFast,nDepth
            print nlines/1000.0,nlines/300.0,fn.find('NFL')
            fid = open(fn.replace('.unp','_cfg.py'),'w')
            fid.write('nSlow=%d\n'%nSlow)
            fid.write('nFast=%d\n'%nFast)
            fid.write('nDepth=%d\n'%nDepth)
            fid.write('nVol=%d\n'%nVol)
            fid.close()
            todelete = fn+'.txt'
            try:
                os.remove(todelete)
            except Exception as e:
                print e,todelete
            
class Processor:

    def __init__(self,fn,dummyServer=False):
        if not fn[-4:]=='.unp':
            sys.exit('Cannot open file %s. It is not a .unp file.'%fn)
            
        # make a working directory:
        self.fulltag = fn.replace('.unp','')
        self.workingDirectory = self.fulltag + '/'
        if not os.path.exists(self.workingDirectory):
            os.makedirs(self.workingDirectory)
            print 'Creating working directory %s'%self.workingDirectory
        else:
            print 'Existing working directory %s will be used'%self.workingDirectory

        # make a log directory:
        self.logDirectory = os.path.join(self.workingDirectory,'log')
        if not os.path.exists(self.logDirectory):
            os.makedirs(self.logDirectory)
            print('Creating directory %s'%self.logDirectory)
        else:
            print('Existing log directory %s will be used'%self.logDirectory)


        self.mlog = Logger(ocfg.MASTER_LOG)
        self.logfn = os.path.join(self.logDirectory,'log.txt')
        self.optlogfn = os.path.join(self.logDirectory,'optimization_log.txt')

        self.log = Logger(self.logfn)
        self.optlog = Logger(self.optlogfn,timestamp=False)


        self.isNfl = fn.lower().find('nfl')>-1
        
        try:
            # find the matching DATA_ROOT
            for dr in ocfg.DATA_ROOTS:
                startidx = self.fulltag.find(dr)
            if startidx>-1:
                branch = self.fulltag[len(dr):]
                if branch.find('\\')>-1:
                    atts = branch.split('\\')
                else:
                    atts = branch.split('/')

            self.attributes = {}
            # atts[0] is the date
            datestr = atts[0]
            if datestr.find('.')>-1:
                datelist = datestr.split('.')
            elif datestr.find('-')>-1:
                datelist = datestr.split('-')
            elif datestr.find('_')>-1:
                datelist = datestr.split('_')
            else:
                datelist = [datestr[:4],datestr[4:6],datestr[6:8]]

            dateintlist = []
            for item in datelist:
                dateintlist.append(int(item))

            # sanity check order of datelist items:
            if dateintlist[0]<dateintlist[1] or dateintlist[0]<dateintlist[2]:
                dateintlist = [0,0,0]
                self.log.log('Trouble parsing date %s. Using 0/0/0.'%datestr)

            self.attributes['year'] = dateintlist[0]
            self.attributes['month'] = dateintlist[1]
            self.attributes['day'] = dateintlist[2]

            subjectName = atts[1]
            subjectName = subjectName.replace(', ','_')
            subjectName = subjectName.replace(',','_')
            subjectName = subjectName.replace(' ','_')
            self.attributes['name'] = subjectName
            self.attributes['eye'] = atts[2].replace(' ','_')

            info = ''
            for item in atts[3:]:
                info = info + item + ' '

            info = info[:-1]
            self.attributes['info'] = info.replace(' ','_')

            
        except Exception as e:
            self.log.log(e)
            self.log.log('Using default attributes.')
            self.attributes['year'] = 1900
            self.attributes['month'] = 1
            self.attributes['day'] = 1
            self.attributes['name'] = 'default_name'
            self.attributes['eye'] = 'NA'
            self.attributes['info'] = 'default_attributes'

            

        self.server = Server(self.log,dummyServer)

        try:
            self.FLIPPED = ocfg.FLIPPED
        except Exception as e:
            print('FLIPPED not found in octtools_cfg; defaulting to False. If b-scans are upside down, please add FLIPPED = True to octtools_cfg.py.')
            self.FLIPPED = False

        ct = ConfigurationTranslator(self.log)
        cfg = ct.getCfg(fn)

        try:
            self.system_id = cfg.system_id
        except Exception as e:
            self.log.log('Cannot get system_id from dataset configuration file. Defaulting to 1G-AO-OCT.')
            self.system_id = '1G-AO-OCT'

        # make a diagnostic png directory:
        self.diagnosticPngDirectory = os.path.join(self.workingDirectory,'diagnostics_pngs')
        if not os.path.exists(self.diagnosticPngDirectory):
            os.makedirs(self.diagnosticPngDirectory)
            self.log.log('Creating directory %s'%self.diagnosticPngDirectory)
        else:
            self.log.log('Existing diagnostic directory %s will be used'%self.diagnosticPngDirectory)

        # make a backup directory:
        self.backupDirectory = os.path.join(self.workingDirectory,'bak')
        if not os.path.exists(self.backupDirectory):
            os.makedirs(self.backupDirectory)
            self.log.log('Creating directory %s'%self.backupDirectory)
        else:
            self.log.log('Existing backup directory %s will be used'%self.backupDirectory)

        self.hdf5fn = fn.replace('.unp','')+'.hdf5'
        try:
            self.h5 = h5py.File(self.hdf5fn,'r+')
        except Exception as e:
            print e
            self.h5 = h5py.File(self.hdf5fn,'w')


        self.nBytes = os.stat(fn).st_size
        
        # these are params for loading the raw data, but it changes shape through
        # the processing steps, and the data size should always be checked on the fly
        self.nSlow = cfg.nSlow
        self.nFast = cfg.nFast
        self.nDepth = cfg.nDepth

        try:
            self.nVol = cfg.nVol
        except Exception as e:
            expectedPixelsPerVolume = self.nSlow*self.nFast*self.nDepth
            self.nVol = self.nBytes/expectedPixelsPerVolume/2 # 2 bytes per pixel

        try:
            self.dispersionCoefs = self.h5['/dispersionCoefs']
            self.log.log('Setting dispersion coefs from %s.'%(self.hdf5fn))
            self.log.log('\t Dispersion coefs: %0.3e and %0.3e.'%(self.dispersionCoefs[0],self.dispersionCoefs[1]))
        except Exception as e:
            print e
            try:
                self.dispersionCoefs = cfg.dispersionCoefs
                self.log.log('Setting dispersion coefs from local config file.')
                self.log.log('\tDispersion coefs: %0.3e and %0.3e.'%(self.dispersionCoefs[0],self.dispersionCoefs[1]))
            except Exception as e:
                try:
                    dispfn = os.path.join(self.workingDirectory,'dispersion.txt')
                    self.dispersionCoefs = np.loadtxt(dispfn)
                    self.log.log('Setting dispersion coefs from %s.'%dispfn)
                    self.log.log('\tDispersion coefs: %0.3e and %0.3e.'%(self.dispersionCoefs[0],self.dispersionCoefs[1]))
                except Exception as e:
                    self.dispersionCoefs = ocfg.DISPERSION_COEFS
                    self.log.log('Setting dispersion coefs from octtools_config.')
                    self.log.log('\tDispersion coefs: %0.3e and %0.3e.'%(self.dispersionCoefs[0],self.dispersionCoefs[1]))

        try:
            self.fastMin = ocfg.FAST_MIN
        except Exception as e:
            self.fastMin = 0

        try:
            self.fastMax = ocfg.FAST_MAX
        except Exception as e:
            self.fastMax = self.nFast


        try:
            self.depthMin = ocfg.DEPTH_MIN
        except Exception as e:
            self.depthMin = 0

        try:
            self.depthMax = ocfg.DEPTH_MAX
        except Exception as e:
            self.depthMax = self.nDepth


        self.path,self.filename = os.path.split(fn)
        self.unpfn = fn
        self.dcfn = self.unpfn.replace('.unp','_dc.npy')
        self.fid = open(fn,'rb')

        try:
            # new method: permit arbitrary order polynomial coefficients, as long
            # as they're ordered by numpy.polyval order
            mappingPolynomial = ocfg.MAPPING_POLYNOMIAL
            self.L = np.polyval(mappingPolynomial,np.arange(self.nDepth))
            self.L_MIN,self.L_MAX = self.L.min(),self.L.max()
        except Exception as e:
            # if that doesn't work, try to load L_MIN and L_MAX explicitly
            print e
            self.L_MIN,self.L_MAX = ocfg.L_MIN,ocfg.L_MAX
            self.L = np.linspace(self.L_MIN,self.L_MAX,self.nDepth)

        

        self.k_in = (2.0*np.pi)/self.L
        self.k_out = np.linspace(self.k_in[0],self.k_in[-1],self.nDepth)
        
        
        
        # self.k = np.linspace(2*np.pi/self.L_MIN,2*np.pi/self.L_MAX,self.nDepth)
        # self.L_out = 2*np.pi/self.k
        # # protect against L_out values outside the L range, due to roundoff:
        # self.L_out[0] = self.L[0]
        # self.L_out[-1] = self.L[-1]


        

        

        vidfn = os.path.join(self.workingDirectory,'vid.txt')
        try:
            fid = open(vidfn)
            self.vid = fid.read()
            fid.close()
            self.log.log('Using unique ID (md5 sum) %s for %s.'%(self.vid,fn))
        except Exception as e:
            print e
            self.log.log('Computing unique ID (md5 sum) for %s.'%fn)
            # make a unique ID based on the md5 sum, in 128MB chunks
            chunkSize = 1024*1024*128
            md5 = hashlib.md5()
            fid = open(fn)
            chunk = fid.read(chunkSize)
            totalBytes = 0
            chunkCount = 0
            while chunk:
                self.log.log('\tchunk %d, total bytes %d'%(chunkCount,totalBytes))
                md5.update(chunk)
                chunk = fid.read(chunkSize)
                chunkCount = chunkCount + 1
                totalBytes = totalBytes + chunkSize
            self.vid = md5.hexdigest()+'_'+self.getSubjectInitials()+'_'+self.getInfo()
            self.log.log('ID: %s'%(self.vid))

            fid.close()
            fid = open(vidfn,'w')
            fid.write(self.vid)
            fid.close()


        self.eid = self.getExperimentId()
        eidfn = os.path.join(self.workingDirectory,'eid.txt')
        fid = open(eidfn,'w')
        fid.write(self.eid)
        fid.close()

        self.sid = self.getSubjectId()
        sidfn = os.path.join(self.workingDirectory,'sid.txt')
        fid = open(sidfn,'w')
        fid.write(self.sid)
        fid.close()

        self.mlog.log('Volume ID: %s'%(self.vid))
        self.mlog.log('Subject ID: %s'%(self.sid))
        self.mlog.log('Experiment ID: %s'%(self.eid))
        
        for key in self.attributes.keys():
            self.mlog.log('%s: %s'%(key,self.attributes[key]))
            
        # now we have to add the following information to the database:
        # volume id, subject id, experiment id, date, eye, info, processed = 0
        self.server.addVolume(self.vid,self.sid,self.eid,self.getDate(),self.getEye(),self.getInfo(),0)

        relputplace = self.sid + '/' + self.eid + '/' + self.vid + '/'
        self.putPlace = os.path.join(ocfg.SERVER_WWW_ROOT, relputplace)
        self.getSubjectInitials()


    
        
    def put(self,src):
        srcpath,srcfn = os.path.split(src)
        dest = os.path.join(self.putPlace,srcfn)
        self.server.put(src,dest,relative=False)
        
    def putLogs(self):
        self.put(self.logfn)
        self.put(self.optlogfn)

    def getDate(self):
        return '%04d-%02d-%02d'%(self.attributes['year'],self.attributes['month'],self.attributes['day'])

    def getEye(self):
        return self.attributes['eye']

    def getInfo(self):
        return self.attributes['info']

    def getSubjectName(self):
        return self.attributes['name']

    def getSubjectInitials(self):
        n = self.attributes['name']
        nlist = n.split('_')
        try:
            initials = nlist[1][:3] + nlist[0][:3]
        except Exception as e:
            print e
            initials = nlist[0][:3]
        return initials
            
    def getExperimentId(self):
        md5 = hashlib.md5()
        md5.update(self.getSubjectName())
        md5.update(self.getDate())
        return md5.hexdigest()

    def getSubjectId(self):
        md5 = hashlib.md5()
        md5.update(self.getSubjectName())
        return md5.hexdigest()


    def deleteDiagnostics(self,backup=True):
        if backup:
            misc.backupDirectory(self.diagnosticPngDirectory)
        misc.deleteDirectoryContents(self.diagnosticPngDirectory)

    def deleteLogs(self,backup=True):
        if backup:
            misc.backupDirectory(self.logDirectory)
        misc.deleteDirectoryContents(self.logDirectory)


    def getFrame(self,volumeIndex,frameIndex,dcSubtract=True,adaptiveDC=False):
        position = volumeIndex * self.nDepth * self.nFast * self.nSlow * 2 + frameIndex * self.nDepth * self.nFast * 2
        self.fid.seek(position,0)
        data = np.fromfile(self.fid,dtype=np.uint16,count=self.nDepth*self.nFast)
        data = data.reshape(self.nFast,self.nDepth).T
        # crop in fast dimension first, to spare processing junk
        data = data[:,self.fastMin:self.fastMax]

        if dcSubtract and adaptiveDC:
            bias = np.convolve(np.ones(data.shape[1]),np.ones(50),'same')
            dc = convolve2d(data,np.ones((1,50)),'same')/bias
            data = data - dc
        elif dcSubtract:
            dc = data.mean(axis=1)
            data = (data.T - dc).T

        return np.flipud(data)


    def saveFrame(self,fn,frame,title=None):
        plt.figure()
        plt.imshow(frame)
        plt.colorbar()
        if not title is None:
            plt.title(title)
        plt.savefig(fn)
        plt.close()

    def depthCentroid(self,frame):
        sy,sx = frame.shape
        xx,yy = np.meshgrid(np.arange(sx),np.arange(sy))
        return np.sum(frame*yy)/np.sum(frame)
        

    def optimizeDispersionCoefficientsDebug(self,N=40,width=10,plot=False):
        return self.getOptimizedDispersionCoefficients(N,width,plot,True)

    def optimizeDispersionCoefficients(self,N=10,width=20,plot=False):
        coefs = self.getOptimizedDispersionCoefficients(N,width,plot,False)
        self.dispersionCoefs = coefs
        self.writeDispersionCoefsToH5(coefs)

    def writeGlobalDispersionCoefs(self,coefs):
        path,fn = os.path.split(self.fulltag)
        dcfn = os.path.join(path,'dispersionCoefs.txt')
        np.savetxt(dcfn,coefs)
        
    def readGlobalDispersionCoefs(self,nfl=False):
        coefs = None
        try:
            path,fn = os.path.split(self.fulltag)
            if nfl:
                dcfn = os.path.join(path,'dispersionCoefsNFL.txt')
            else:
                dcfn = os.path.join(path,'dispersionCoefs.txt')
            self.log.log('\t...%s'%dcfn)
            coefs = np.loadtxt(dcfn)
        except Exception as e:
            print e, 'returning None'
        return coefs

    def readDispersionCoefsFromH5(self):
        try:
            return self.h5['dispersionCoefs']
        except Exception as e:
            raise e

    def writeDispersionCoefsToH5(self,coefs):
        try:
            self.log.log('Trying to write dispersion coefficients to H5 file.')
            self.h5['dispersionCoefs'] = coefs
            self.log.log('\tSuccess (1).')
        except Exception as e:
            print e
            try:
                self.log.log('Trying to write dispersion coefficients to H5 file.')
                self.h5['dispersionCoefs'][...] = coefs
                self.log.log('\tSuccess (2).')
            except Exception as ee:
                print ee
                try:
                    self.log.log('Trying to write dispersion coefficients to H5 file.')
                    del self.h5['dispersionCoefs']
                    junk = self.h5.create_dataset('dispersionCoefs',data=coefs)
                    self.log.log('\tSuccess (3).')
                except Exception as eee:
                    sys.exit(eee)


    def makeTestFrame(self,N,width):

        temp = self.getFrame(0,0)
        sy,sx = temp.shape
        testFrame = np.zeros((sy,N*width))
        for n in range(N):
            vidx = np.random.randint(self.nVol)
            sidx = np.random.randint(self.nSlow)
            temp = self.getFrame(vidx,sidx)
            temp = temp[:,sx/2:sx/2+width]
            testFrame[:,n*width:(n+1)*width] = temp
        
           
        return testFrame


    def showDispersionCompensation(self,coefs=None,showPlot=False):
        if coefs is None:
            coefs = self.dispersionCeofs
        #tf = self.makeTestFrame(100,10)
        tf = self.getFrame(0,round(self.nSlow/3))
        pre = np.abs(self.processFrame(tf,[0.0,0.0,0.0,0.0]))
        post = np.abs(self.processFrame(tf,coefs))
        presnr = pre.max()/self.noiseRMS(pre)
        postsnr = post.max()/self.noiseRMS(post)

        plt.figure(figsize=(8,4))
        plt.subplot(1,2,1)
        plt.imshow((np.abs(pre)))
        plt.title('pre compensation SNR = %0.1f'%presnr)
        plt.colorbar()
        plt.subplot(1,2,2)
        plt.imshow((np.abs(post)))
        plt.title('post compensation SNR = %0.1f'%postsnr)
        plt.colorbar()
        cfn = os.path.join(self.diagnosticPngDirectory,'dispersionComparison.png')
        plt.savefig(cfn)
        self.put(cfn)

        if showPlot:
            pass
            #plt.show()

        plt.close()

    def getOptimizedDispersionCoefficients(self,N,width,plot,debug):
        testFrame = self.makeTestFrame(N,width)
        
        obj = lambda c: self.dispersionObjective(testFrame,c,plot=plot)

#        self.log.log('Getting initial values for optimization from...')
#        c0 = self.readGlobalDispersionCoefs(self.isNfl)
#        if c0 is None:
#            c0 = self.dispersionCoefs[:2]

        c0 = [0.0,0.0]

        precoef = [0.,0.,0.,0.]
        postcoef = [0.,0.,0.,0.]

        #lowers = ocfg.DISPERSION_COEFS_LBOUNDS[:2]
        #uppers = ocfg.DISPERSION_COEFS_UBOUNDS[:2]
        lowers = [c0[0]-2e-16,c0[1]-2e-10]
        uppers = [c0[0]+2e-16,c0[1]+2e-10]
        
        bounds3 = (lowers[0],uppers[0])
        bounds2 = (lowers[1],uppers[1])


        if plot:
            plt.figure()

        t0 = time()
        self.log.log('optimization start')
        self.optlog.log('# optimization start')
        self.optlog.log('# starting with coefs %0.3e,%0.3e'%(c0[0],c0[1]))
        self.optlog.log('# order3\torder2\t\tmetric')


        method='brute'

        if method=='brute':
            result = optimize.brute(obj,(bounds3,bounds2),Ns=41,finish=None)
            bounds3a = (result[0]-1e-17,result[0]+1e-17)
            bounds2a = (result[1]-1e-11,result[1]+1e-11)
            result = optimize.brute(obj,(bounds3a,bounds2a),Ns=11,finish=None)
        elif method=='anneal':
            result = optimize.anneal(obj,x0=c0,schedule='fast',lower=lowers,upper=uppers,maxeval=None)
        else:
            sys.exit('Please provide an optimization method.')
            
        t_elapsed = time() - t0

        if plot:
            plt.ioff()
            plt.close()

        if method=='anneal':
            posttemp = result[0]
        elif method=='brute':
            posttemp = result
        else:
            sys.exit('Please provide an optimization method.')

        postcoef[0] = posttemp[0]
        postcoef[1] = posttemp[1]

        self.plotOptLog()
        self.log.log('optimization time elapsed: %f'%t_elapsed)
        self.optlog.log('# optimization time elapsed: %f'%t_elapsed)
        self.showDispersionCompensation(coefs=postcoef,showPlot=False)
        return postcoef

    def computeIntensitySquared(self,frame):
        frame = frame**2.0
        maxes = frame.max(axis=0)
        return maxes.mean()

    def computeNonsense(self,frame):
        #fnoise = np.std(frame[-50:,:],axis=0)
        fmax = np.max(frame,axis=0)
        fmaxidx = np.argsort(fmax)
        SNR = np.median(fmax[fmaxidx[:-50]])
        return SNR

    def computeSNR(self,frame,mode='maximum'):
        frame = frame
        vprof = frame.mean(axis=1)
        noiseidx = np.argsort(vprof)[:20]
        noiserms = frame[noiseidx,:].std()
        fmax = np.max(frame,axis=0)
        if mode=='maximum':
            SNR = np.max(fmax/noiserms)
        elif mode=='median':
            SNR = np.median(fmax/noiserms)
        else:
            sys.exit('computeSNR: invalid mode')

        return SNR

    def computeContrastMax(self,frame):
        pmax = np.max(frame,axis=0)
        pmin = np.median(frame,axis=0)
        contrast = (pmax-pmin)/(pmax+pmin)
        return contrast.max()*100.0

    def computeContrast(self,frame):
        pmax = np.max(frame,axis=0)
        pmin = np.median(frame,axis=0)
        contrast = (pmax-pmin)/(pmax+pmin)
        return contrast.mean()*100.0
        
    def computeProfileContrast(self,frame):
        frame = frame
        p = frame.mean(axis=1)
        return (p.max()-p.min())/(p.max()+p.min())*100.0
        
    def computeMaxGradients(self,frame):
        return np.std(np.diff(frame,axis=0),axis=0).max()
        
    def computeHighPix(self,frame):
        #frame = frame
        fmean = frame.mean()
        fstd = frame.std()
        thresh = fmean+3*fstd
        #thresh = 15000
        sy,sx = frame.shape
        return len(np.where(frame>thresh)[0])/float(sx)/float(sy)
    
    def noiseRMS(self,frame):
        if np.iscomplex(frame[0,0]):
            frame = np.abs(frame)
        sy,sx = frame.shape
        prof = np.mean(frame,axis=1)
        nnoise = int(0.1*sy)
        noiseIdx = np.argsort(prof)[:nnoise]
        noise = frame[noiseIdx,:]
        return noise.std()


    def dispersionObjective(self,frame,coefs,plot=False):
        allcoefs = [coefs[0],coefs[1],0.0,0.0]
        frame = np.abs(self.processFrame(frame,allcoefs,plot=plot))
        frame = frame[:-50,:]
        Isquare = self.computeIntensitySquared(frame)
        self.optlog.log('%0.3e\t%0.3e\t%0.4e'%(coefs[0],coefs[1],Isquare))
        return 1.0/Isquare

    def dispersionObjective1(self,frame,coefs,plot=False):
        allcoefs = [coefs[0],coefs[1],0.0,0.0]
        frame = np.abs(self.processFrame(frame,allcoefs,plot=plot))
        frame = frame[:-50,:]
        SNR = self.computeSNR(frame)
        self.optlog.log('%0.3e\t%0.3e\t%0.3e'%(coefs[0],coefs[1],SNR))
        return 1.0/SNR

    def dispersionObjective0(self,frame,coefs,plot=False):
        allcoefs = [coefs[0],coefs[1],0.0,0.0]
        frame = np.abs(self.processFrame(frame,allcoefs,plot=plot))
        c = self.computeContrast(frame)
        self.optlog.log('%0.3e\t%0.3e\t%0.4f'%(coefs[0],coefs[1],c))
        return 1.0/c


    def plotOptLog(self,showPlot=False):
        olog = np.loadtxt(self.optlogfn)
        c3 = olog[:,0]
        c2 = olog[:,1]
        snr = olog[:,2]

        fig = plt.figure()
        ax = fig.gca(projection='3d')

        X = c3#np.arange(-5, 5, 0.25)
        Y = c2#np.arange(-5, 5, 0.25)
        Z = snr
        surf = ax.scatter(X, Y, Z)
        ax.set_xlabel('3rd order')
        ax.set_ylabel('2nd order')
        ax.set_zlabel('SNR ($I^2$)')
        outfn = os.path.join(self.diagnosticPngDirectory,'dispersionSurface.png')
        plt.title('dispersion surface')
        plt.savefig(outfn)
        self.put(outfn)

        if showPlot:
            pass
            #plt.show()        

        plt.close()


        plt.figure()
        t = np.arange(len(c3))
        ax1 = plt.subplot(2,1,1)
        ax1.plot(t,c3,'k.')
        ax1.set_ylabel('$a$ in $e^{ax^3 + bx^2}$')

        ax2 = plt.twinx()
        ax2.plot(t,snr,'r-')
        ax2.set_ylabel('SNR $I^2$', color='r')
        for tl in ax2.get_yticklabels():
            tl.set_color('r')



        ax3 = plt.subplot(2,1,2)
        ax3.plot(t,c2,'k.')
        ax3.set_ylabel('$b$ in $e^{ax^3 + bx^2}$')


        ax4 = plt.twinx()
        ax4.plot(t,snr,'r-')
        ax4.set_ylabel('SNR $I^2$', color='r')
        for tl in ax4.get_yticklabels():
            tl.set_color('r')

        plt.title('dispersion coefficient optimization')
        outfn = os.path.join(self.diagnosticPngDirectory,'dispersionOptimizationTrace.png')
        plt.savefig(outfn)
        self.put(outfn)

        if showPlot:
            pass
            #plt.show()
        
        plt.close()

    def processFrame(self,frame,dispersionCoefs=None,plot=False):
        # If plot is True, the caller must initialize a figure and set plt.ion().

        if dispersionCoefs is None:
            dispersionCoefs = self.dispersionCoefs

        # map into k-space (right now 1G and 2G systems have opposite pixel-k relationships (i.e. increasing k w/ pixel vs. decreasing k with pixel)
        
        # 2G-AO-OCT working code:
        if self.system_id=='2G-AO-OCT':
            kInterpolator = interpolate.interp1d(self.k_in,frame,axis=0,copy=False)
            frame = kInterpolator(self.k_out)
            
        # 1G-AO-OCT working code:
        if self.system_id=='1G-AO-OCT':
            kInterpolator = interpolate.interp1d(self.k_in[::-1],frame,axis=0,copy=False)
            frame = kInterpolator(self.k_out[::-1])

        #frame = hilbert(frame,axis=0)

        # dispersion compensate
        dispersionAxis = self.k_out - self.k_out.mean()
        phase = np.exp(1j*np.polyval(dispersionCoefs,dispersionAxis))

        frame = frame * phase[None].T
        preframe = np.abs(frame)

        # fourier transform
        frame = np.fft.fftshift(np.fft.fft(frame,axis=0),axes=0)

        frame = frame[self.depthMin:self.depthMax,:]

        if self.FLIPPED:
            frame = np.flipud(frame)
        
        if plot:
            fsnr = self.computeSNR(frame)
            plt.subplot(1,2,1)
            plt.cla()
            plt.imshow(preframe)
            plt.subplot(1,2,2)
            plt.cla()
            plt.imshow(np.abs(frame))
            plt.title(fsnr)
            plt.pause(.1)

        return frame


    def z_align_rough(self,iVol=0,do_plot=True):
        # this method returns a coarsely aligned volume
        # if '/bscan_alignment_vector' exists, it uses this to
        # align the b-scans; if not, it creates that vector first.
        try:
            zvec = self.h5['/bscan_alignment_vector']
        except Exception as e:
            try:
                del self.h5['/bscan_alignment_vector']
            except Exception as e:
                self.log.log('No bscan_alignment_vector to delete from h5. This is good.')
            self.log.log('z_align_1: no bscan_alignment_vector in hdf5 file. Generating.')
            proj = np.mean(np.abs(self.h5['processed'][iVol,:,:,:]),axis=1)
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
            self.h5['/bscan_alignment_vector'] = zvec

        nv,ns,nf,nd = self.h5['processed'].shape
        slow1 = iVol*self.nSlow
        slow2 = slow1 + self.nSlow
        zoff = np.max(zvec)
        
        outvol = np.zeros((ns,nf,nd+zoff),dtype=np.complex64)
        for k in range(ns):
            outvol[k,:,zvec[k]:zvec[k]+nd] = self.h5['processed'][iVol,k,:,:]


        try:
            zvec = self.h5['/fast_axis_alignment_vector']
        except Exception as e:
            try:
                del self.h5['/fast_axis_alignment_vector']
            except Exception as e:
                pass
            self.log.log('No fast axis alignment vector in h5. Generating.')

            proj = np.mean(np.abs(outvol),axis=0)
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
            self.h5['/fast_axis_alignment_vector'] = zvec


        ns,nf,nd = outvol.shape
        zoff = np.max(zvec)
        
        outvol2 = np.zeros((ns,nf,nd+zoff),dtype=np.complex64)
        for k in range(nf):
            outvol2[:,k,zvec[k]:zvec[k]+nd] = outvol[:,k,:]

        return outvol2
            

    def archiveLRPModel(self,model):
        model.write_to_h5(self.h5)
    
    def align(self,plot=False):
        # This method does sequential dewarping of processed b-scans
        # the dewarping is 'virtual', in the sense that the only
        # data stored to the HDF5 file are x,y,z offsets for each A-line.
        # Each offset dimension is stored in a nVol x nSlow x nFast
        # matrix.
        
        try:
            del self.h5['/aligned']
        except Exception as e:
            print e
        
        nDepth = self.depthMax - self.depthMin
        nFast = self.fastMax - self.fastMin
        nSlow = self.nSlow
        nVol = self.nVol

        slowOffset = np.zeros((nVol,nSlow,nFast))
        fastOffset = np.zeros((nVol,nSlow,nFast))
        depthOffset = np.zeros((nVol,nSlow,nFast))
        
    

    def alignTwoBscans(self,ref,tar,plot=False):
        # This method returns two sets of numbers, each equal
        # in size to the width of the TAR b-scan, depthOffsets
        # and fastOffsets. The values in each indicate how many
        # pixels down and to the right to move each a-line in TAR,
        # in order to best aling with REF.
        ref = np.abs(ref)
        tar = np.abs(tar)

        gyshift,gxshift,gpeakVal,gnxcval = self.nxcorr2same(ref,tar)


        if gpeakVal>.3:
            self.log.log('Global inter-frame shifts %d (x), %d (y) with correlation %0.2f.'%(gxshift,gyshift,gpeakVal))
        else:
            self.log.log('Correlation %0.2f < 0.3. Global inter-frame shifts set to 0, 0.'%gpeakVal)
            gyshift = 0
            gxshift = 0

        if gyshift>10 or gxshift>10:
            self.log.log('Large global shift. Verify correlation is reasonable.')
            
        sy,sx = ref.shape
        extRef = np.ones((sy+np.abs(gyshift),sx+np.abs(gxshift)))
        extTar = np.ones((sy+np.abs(gyshift),sx+np.abs(gxshift)))
        

        if gyshift>=0 and gxshift>=0:
            extRef[gyshift:gyshift+sy,gxshift:gxshift+sx] = ref
            extTar[:sy,:sx] = tar
            tarXOff = gxshift
        elif gyshift<0 and gxshift<0:
            extTar[-gyshift:-gyshift+sy,-gxshift:-gxshift+sx] = tar
            extRef[:sy,:sx] = ref
            tarXOff = 0
        elif gyshift>=0 and gxshift<0:
            extRef[gyshift:gyshift+sy,:sx] = ref
            extTar[:sy,-gxshift:-gxshift+sx] = tar
            tarXOff = 0
        else:
            extTar[-gyshift:-gyshift+sy,:sx] = tar
            extRef[:sy,gxshift:gxshift+sx] = ref
            tarXOff = gxshift


        if plot:
            plt.figure(figsize=(12,6))
            plt.ion()

            plt.subplot(2,2,1)
            plt.imshow(extRef)
            plt.draw()

        tsy,tsx = extTar.shape
        rsy,rsx = extRef.shape

        stripWidth = 21

        centers = []
        yshifts = []
        xshifts = []
        peakVals = []


        stepSize = 20.0
        sigma = stepSize*2
        for center in np.arange(0,tsx,stepSize):
            xvec = np.arange(tsx)-center
            g = np.exp(-xvec**2/(2*sigma**2))
            g = g/g.max()
            temp = extTar*g

            #temp = np.zeros(extTar.shape)
            #x1 = center-sigma/2
            #x2 = center+sigma/2
            #temp[:,x1:x2] = extTar[:,x1:x2]

            yshift,xshift,peakVal,nxcim = self.nxcorr2same(extRef,temp,ymax=5,xmax=5) 
            yshifts.append(yshift)
            xshifts.append(xshift)
            peakVals.append(peakVal)
            centers.append(center)
            self.log.log('Windowed shift at %d: %d (x), %d (y), with correlation of %0.2f'%(center,xshift,yshift,peakVal))

            if plot:
                plt.subplot(2,2,2)
                plt.cla()
                plt.imshow(temp)
                plt.subplot(2,2,3)
                plt.cla()
                plt.imshow(nxcim)
                plt.subplot(2,2,4)
                plt.cla()
                maxprof = nxcim.max(axis=0)
                plt.plot(maxprof)
                maxidx = np.argmax(maxprof)
                plt.plot(maxidx,maxprof[maxidx],'ks')
                plt.draw()


        xshifts = np.array(xshifts)
        yshifts = np.array(yshifts)
        peakVals = np.array(peakVals)
        centers = np.array(centers)

        # remove outliers from xshifts and yshifts
        nStd = 10
        xut = np.mean(xshifts)+np.std(xshifts)*nStd
        xlt = np.mean(xshifts)-np.std(xshifts)*nStd
        yut = np.mean(yshifts)+np.std(yshifts)*nStd
        ylt = np.mean(yshifts)-np.std(yshifts)*nStd
        

        xgood = np.logical_and(xshifts<=xut,xshifts>=xlt)
        ygood = np.logical_and(yshifts<=yut,yshifts>=ylt)

        xbad = np.logical_or(xshifts>xut,xshifts<xlt)
        ybad = np.logical_or(yshifts>yut,yshifts<ylt)

        good = np.where(np.logical_and(xgood,ygood))[0]
        bad = np.where(np.logical_or(xbad,ybad))[0]

        if len(bad):
            self.log.log('Shifts at the following positions were not used:')
            self.log.log(bad)

        xshifts = xshifts[good]
        yshifts = yshifts[good]

        if False:
            print good
            print bad

            print xshifts
            print np.logical_and(xshifts<xut,xshifts>xlt)
            print xlt,xut
            print xgood

            print yshifts
            print np.logical_and(yshifts<yut,yshifts>ylt)
            print ylt,yut
            print ygood


        xshifts = xshifts+gxshift
        yshifts = yshifts+gyshift
        useIdx = np.where(peakVals>0.1)[0]

        usePolynomials = False
        xin = np.arange(tsx)
        if usePolynomials:
            px = np.polyfit(centers[useIdx],xshifts[useIdx],3)
            py = np.polyfit(centers[useIdx],yshifts[useIdx],3)
            xfit = np.polyval(px,xin)
            yfit = np.polyval(py,xin)
        else:
            # use interpolation
            xinterp = interpolate.interp1d(centers[useIdx],xshifts[useIdx],'linear',bounds_error=False,fill_value=0)
            yinterp = interpolate.interp1d(centers[useIdx],yshifts[useIdx],'linear',bounds_error=False,fill_value=0)
            xfit = xinterp(xin)
            yfit = yinterp(xin)

        xfit = xfit[tarXOff:tarXOff+sx]
        yfit = yfit[tarXOff:tarXOff+sx]

        if plot:
            plt.close()
            plt.ioff()
        
        return yfit,xfit,gyshift,gxshift


    def clearArchive(self,key):
        try:
            del self.h5[key]
        except Exception as e:
            print e


    def archivePairwiseShifts(self,plot=False):

        try:
            nVol,nSlow,nFast,nDepth = self.h5['/processed'].shape
        except Exception as e:
            self.log.log('Cannot access /processed in %s.'%self.hdf5fn)
            sys.exit(e)

        for key in ['globalXShifts','globalYShifts','localXShifts','localYShifts']:
            self.clearArchive(key)

        globalXShifts = np.zeros((nVol,nSlow))
        globalYShifts = np.zeros((nVol,nSlow))
        localXShifts = np.zeros((nVol,nSlow,nFast))
        localYShifts = np.zeros((nVol,nSlow,nFast))

        for iVol in range(nVol):
            for iSlow in range(1,nSlow):
                self.log.log('Computing pairwise shifts between frames %d and %d.'%(iSlow,iSlow-1))
                ref = self.h5['/processed'][iVol][iSlow-1].T
                tar = self.h5['/processed'][iVol][iSlow].T
                lys,lxs,gys,gxs = self.alignTwoBscans(ref,tar,plot=plot)
                localYShifts[iVol][iSlow] = lys
                localXShifts[iVol][iSlow] = lxs
                globalYShifts[iVol][iSlow] = gys
                globalXShifts[iVol][iSlow] = gxs

                print globalXShifts.mean()
                print localXShifts.mean()
                sleep(2)
            
        self.h5['globalXShifts'] = globalXShifts
        self.h5['globalYShifts'] = globalYShifts
        self.h5['localXShifts'] = localXShifts
        self.h5['localYShifts'] = localYShifts

        np.savetxt('lystemp.txt',localYShifts)
        np.savetxt('lxstemp.txt',localXShifts)
        np.savetxt('gystemp.txt',globalYShifts)
        np.savetxt('gxstemp.txt',globalXShifts)

        


    def interpolateFrame(self,frame,yin,yout,xin,xout):
        #f = interpolate.interp2d(xin,yin,frame)
        return self.bilinear_interpolate(frame,xout,yout)
        #return f(xout,yout)


    def bilinear_interpolate(self, im, x, y):
        x = np.asarray(x)
        y = np.asarray(y)

        x0 = np.floor(x).astype(int)
        x1 = x0 + 1
        y0 = np.floor(y).astype(int)
        y1 = y0 + 1

        x0 = np.clip(x0, 0, im.shape[1]-1);
        x1 = np.clip(x1, 0, im.shape[1]-1);
        y0 = np.clip(y0, 0, im.shape[0]-1);
        y1 = np.clip(y1, 0, im.shape[0]-1);

        Ia = im[ y0, x0 ]
        Ib = im[ y1, x0 ]
        Ic = im[ y0, x1 ]
        Id = im[ y1, x1 ]

        wa = (x1-x) * (y1-y)
        wb = (x1-x) * (y-y0)
        wc = (x-x0) * (y1-y)
        wd = (x-x0) * (y-y0)

        return wa*Ia + wb*Ib + wc*Ic + wd*Id


    def archiveProcessed(self,plot=False):
        if plot:
            plt.figure()

        nFast = self.fastMax - self.fastMin
        nDepth = self.depthMax - self.depthMin
        nSlow = self.nSlow
        nVol = self.nVol

        pmax = -np.inf
        pmin = np.inf

        try:
            del self.h5['/processed']
        except Exception as e:
            print e

        self.log.log('Creating /processed matrix in hdf5 file %s. Please wait.'%(self.hdf5fn))
        processed = self.h5.create_dataset('/processed',(nVol,nSlow,nFast,nDepth),dtype=np.complex64)

        for iVol in range(nVol):
            for iSlow in range(nSlow):
                self.log.log('Processing volume %02d of %02d, frame %03d of %03d.'%(iVol+1,nVol,iSlow+1,nSlow))
                frame = self.getFrame(iVol,iSlow)
                frame = self.processFrame(frame,plot=plot)
                frame = frame.T
                if plot:
                    plt.imshow(np.abs(frame))
                    plt.pause(.1)
                processed[iVol,iSlow,:,:] = frame
                if np.abs(frame).max()>pmax:
                    pmax = np.abs(frame).max()
                if np.abs(frame).min()<pmin:
                    pmin = np.abs(frame).min()

        try:
            self.h5['/processed_max'] = pmax
        except Exception as e:
            self.log.log('Cannot write processed_max to h5 file. Trying to delete and rewrite.')
            del self.h5['/processed_max']
            self.h5['/processed_max'] = pmax
            
        try:
            self.h5['/processed_min'] = pmin
        except Exception as e:
            self.log.log('Cannot write processed_min to h5 file. Trying to delete and rewrite.')
            del self.h5['/processed_min']
            self.h5['/processed_min'] = pmin

        if plot:
            plt.close()

    def makeMovie(self):

        try:
            nVol,nSlow,nFast,nDepth = self.h5['processed'].shape
        except Exception as e:
            print('Processor.makeMovie: %s'%e)
            return

        testend = min(10,nSlow)
        testchunk = np.abs(self.h5['processed'][0,0:testend:2,:,:])
        cmin = np.median(testchunk)  - testchunk.std()
        cmax = cmin + 7*testchunk.std()

        avifn = os.path.join(self.workingDirectory,'processed.avi')
        mov = Movie(avifn,cmin = cmin, cmax = cmax, make_wmv=False)

        for iVol in range(nVol):
            for iSlow in range(nSlow):
                print('Adding frame %d of %d.'%(iVol*nSlow+iSlow,nVol*nSlow))
                frame = self.h5['processed'][iVol,iSlow,:,:].T
                frame = np.abs(frame)
                #imh = plt.imshow(frame)
                #imh.set_clim([cmin,cmax]) 
                #plt.draw()
                mov.add(frame)

        mov.make()
        self.put(avifn.replace('.avi','.webm'))
        self.put(avifn.replace('.avi','_preview.png'))
        self.put(avifn.replace('.avi','_thumb.png'))

class Logger:

    def __init__(self,fn='./log.txt',overwrite=False,backup=True,timestamp=True):
        self.timestamp = timestamp
        self._fn = fn
        self.log('# Initializing logging to %s.'%self._fn)
        if overwrite:
            nowStr = datetime.datetime.now().strftime("%Y.%m.%d.%H.%M.%S")
            if os.path.exists(fn):
                if backup:
                    shutil.move(fn,fn.replace('.txt','.bak.%s'%nowStr))
                else:
                    os.remove(fn)
        
    def log(self,text,newLine=True):
        if newLine:
            if self.timestamp:
                nowStr = misc.nowStr()#datetime.datetime.now().strftime("%Y.%m.%d.%H.%M.%S")
                outstr = '%s\t%s'%(nowStr,text)
            else:
                outstr = '%s'%text
            print outstr
        else:
            outstr = '%s'%text

        nl = '%s'%'\n'
        self._fid = open(self._fn,'a')
        self._fid.write(outstr)
        if newLine:
            self._fid.write(nl)
        else:
            self._fid.write(' ')
        self._fid.close()



