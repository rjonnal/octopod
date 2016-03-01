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
