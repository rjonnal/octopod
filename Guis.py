import sys,os
from PySide import QtGui, QtCore
import OCTTools
import numpy as np
from matplotlib import pyplot as plt
import octtools_config


class PictureLabel(QtGui.QLabel):

    pictureClicked = QtCore.Signal(list) # can be other types (list, dict, object...)
    pictureDoubleClicked = QtCore.Signal(list) # can be other types (list, dict, object...)

    def __init__(self, parent=None):
        super(PictureLabel, self).__init__(parent)        

    def mouseDoubleClickEvent(self, event):
        print "from PictureLabel.mouseDoubleClickEvent"
        self.pictureDoubleClicked.emit(event.pos())

    def mousePressEvent(self, event):
        print "from PictureLabel.mousePressEvent"
        width = float(self.size().width())
        height = float(self.size().height())
        x = event.pos().x()
        y = event.pos().y()
        self.pictureClicked.emit(event.pos())

class ProcessingWindow(QtGui.QWidget):

    def __init__(self, processor=None):

        super(ProcessingWindow, self).__init__()

        # openFile = QtGui.QAction(QtGui.QIcon('open.png'), 'Open', self)
        # openFile.setShortcut('Ctrl+O')
        # openFile.setStatusTip('Open new processor')
        # openFile.triggered.connect(self.showDialog)

        # menubar = self.menuBar()
        # fileMenu = menubar.addMenu('&File')
        # fileMenu.addAction(openFile)   


        #print processor
        if processor is None:
            processor = self.showDialog()


        self.initializeProcessor(processor)

        self._layout = QtGui.QGridLayout()


        #self.rcanvas = QtGui.QLabel(self)
        #self.rcanvas.setText('No image.')
        #self.rcanvas.setMinimumWidth(600)

        self.pcanvas = PictureLabel(self)
        self.pcanvas.setText('No image.')
        self.pcanvas.setMinimumWidth(600)

        self.pcanvas.pictureDoubleClicked.connect(self.resetImage)
        self.pcanvas.pictureClicked.connect(self.zoomImage)


        #self._layout.addWidget(self.rcanvas, 0, 0)
        self._layout.addWidget(self.pcanvas, 0, 0)

        #splitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
        #splitter.addWidget(self.spinBox('low',self.setRawLow,self.getRawLow,1))
        #splitter.addWidget(self.spinBox('high',self.setRawHigh,self.getRawHigh,1))
        #self._layout.addWidget(splitter, 1, 0)


        vsplitter = QtGui.QSplitter(QtCore.Qt.Vertical)
        splitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
        splitter.addWidget(self.spinBox('low',self.setProcLow,self.getProcLow,100))
        splitter.addWidget(self.spinBox('high',self.setProcHigh,self.getProcHigh,100))
        splitter.addWidget(self.checkBox('log',self.toggleLogMode,self.getLogMode))


        self.infoLabel = QtGui.QLabel(self)
        self.infoLabel.setText('No info.')

        vsplitter.addWidget(self.infoLabel)
        vsplitter.addWidget(splitter)

        self._layout.addWidget(vsplitter, 1, 0)


        self.c3multiplier = 1e-18
        self.c3slider = QtGui.QSlider()
        self.c3slider.setRange(-1000,1000)
        self.c3slider.setValue(self.dispersionCoefs[0]/self.c3multiplier)
        self.c3slider.setTickPosition(QtGui.QSlider.TicksLeft)
        self.c3slider.setTickInterval(100)
        self.c3slider.setTracking(False)
        self.c3slider.valueChanged[int].connect(self.setC3)

        self.c2multiplier = 1e-12
        self.c2slider = QtGui.QSlider()
        self.c2slider.setRange(-1000,1000)
        self.c2slider.setValue(self.dispersionCoefs[1]/self.c2multiplier)
        self.c2slider.setTickPosition(QtGui.QSlider.TicksLeft)
        self.c2slider.setTickInterval(100)
        self.c2slider.setTracking(False)
        self.c2slider.valueChanged[int].connect(self.setC2)



        #splitter.addWidget(self.c3slider)
        #self._layout.addWidget(splitter,0,2)

        self._layout.addWidget(self.c3slider,0,2)
        self._layout.addWidget(self.c2slider,0,3)


        splitter = QtGui.QSplitter(QtCore.Qt.Vertical)
        splitter.addWidget(self.button('Open',self.openDataSet))
        splitter.addWidget(self.button('Update',self.regetAndUpdate))
        splitter.addWidget(self.button('Write coefs',self.writeCoefs))
        splitter.addWidget(self.button('Write global coefs',self.writeGlobalCoefs))
        splitter.addWidget(self.button('Set to global',self.setToGlobal))
        splitter.addWidget(self.button('Print',self.printImage))
        splitter.addWidget(self.spinBox('iVol',self.setIVol,self.getIVol,1,aMax=self.processor.nVol-1,aMin=0))
        splitter.addWidget(self.spinBox('iSlow',self.setISlow,self.getISlow,1,aMax=self.processor.nSlow-1,aMin=0))

        self._layout.addWidget(splitter,0,6)


        self._layout.setRowStretch(0,1)
        self._layout.setRowStretch(1,0)
        self._layout.setRowStretch(2,0)
        self._layout.setRowMinimumHeight(0,400)
        self.setLayout(self._layout)

        self.imageStatsStringBuffer = ['']*5

        test = self.getRawImage(True)
        self.sy,self.sx = test.shape
        self.rqi = QtGui.QImage(test.data,self.sx,self.sy,self.sx,QtGui.QImage.Format_Indexed8)

        test = self.getProcessedImage(True)
        self.sy,self.sx = test.shape
        self.pqi = QtGui.QImage(test.data,self.sx,self.sy,self.sx,QtGui.QImage.Format_Indexed8)

        
        self.aspect = QtCore.Qt.KeepAspectRatio
        self.setGeometry(30,30,1000,800)
        self.scaleFactor = 1.0
        self.infoLabel.setText(self.getInfoString())
        self.path = octtools_config.DATA_ROOTS[0]
        self.show()
        self.update(True)

    def blahblah(self,val):
        print 'blah blah %d'%val
        
    def getImageStatsString(self,frame):
        sy,sx = frame.shape
        imstd = frame.std()
        immean = frame.mean()
        immedian = np.median(frame)
        immax = frame.max()
        immin = frame.min()
        imsum = np.sum(frame)
        
        gradient = np.sqrt(np.sum(np.diff(frame,axis=0)**2)/float(sy))
        outstr = 'M=%0.2e\tm=%0.2e\tu=%0.2e\td=%0.2e\tu-d=%0.2e\tS=%0.2e\tg=%0.2e'%(immax,immin,immean,immedian,immean-immedian,imsum,gradient)
        return outstr
        
    def printImage(self):
        plt.figure(figsize=(7,4))
        ih = plt.imshow(self.proc)
        ih.set_clim((self.procLow,self.procHigh))
        plt.colorbar()
        plt.show()
        
    def zoomImage(self):
        self.scaleFactor *= 1.5
        #self.pcanvas.resize(self.scaleFactor * self.pcanvas.pixmap().size())

        #self.adjustScrollBar(self.scrollArea.horizontalScrollBar(), factor)
        #self.adjustScrollBar(self.scrollArea.verticalScrollBar(), factor)

        #self.zoomInAct.setEnabled(self.scaleFactor < 3.0)
        #self.zoomOutAct.setEnabled(self.scaleFactor > 0.333)

    def resetImage(self):
        self.scaleFactor = 1
        #self.pcanvas.resize(self.scaleFactor * self.pcanvas.pixmap().size())

    def showDialog(self,path=None):
        if path is None:
            path = octtools_config.DATA_ROOTS[0]
        
        fname, _ = QtGui.QFileDialog.getOpenFileName(self, 'Open file',
                    path,'*.unp')
        path,f = os.path.split(fname)
        self.path = path
        return OCTTools.Processor(fname,True)
        
    def openDataSet(self):
        p = self.showDialog(self.path)
        self.initializeProcessor(p)
        test = self.getRawImage(True)
        self.sy,self.sx = test.shape
        self.rqi = QtGui.QImage(test.data,self.sx,self.sy,self.sx,QtGui.QImage.Format_Indexed8)

        test = self.getProcessedImage(True)
        self.sy,self.sx = test.shape
        self.pqi = QtGui.QImage(test.data,self.sx,self.sy,self.sx,QtGui.QImage.Format_Indexed8)

        self.infoLabel.setText(self.getInfoString())
        self.show()
        self.update(True)

    def initializeProcessor(self,processor):
        # contrast values
        self.rawHigh = 200.0
        self.rawLow = -200.0
        self.procHigh = 75000.0
        self.procLow = 7500.0

        self.limitMultiplier = 10.0

        self.logMode = False

        self.processor = processor

        self.iVol = 0
        self.iSlow = 0
        self.dispersionCoefs = [0.0,0.0,0.0,0.0]
        try:
            self.dispersionCoefs[:] = self.processor.readDispersionCoefsFromH5()
        except Exception as e:
            print e


    def setIVol(self,n):
        self.iVol = n
        self.update(True)

    def setISlow(self,n):
        self.iSlow = n
        self.update(True)

    def getISlow(self):
        return self.iSlow

    def getIVol(self):
        return self.iVol


    def writeCoefs(self):
        self.processor.writeDispersionCoefsToH5(self.dispersionCoefs)

    def writeGlobalCoefs(self):
        self.processor.writeGlobalDispersionCoefs(self.dispersionCoefs)

    def readGlobalCoefs(self):
        return self.processor.readGlobalDispersionCoefs()
        
    def setToGlobal(self):
        c = self.readGlobalCoefs()
        if c is not None:
            self.dispersionCoefs = c
            self.c3slider.setValue(self.dispersionCoefs[0]/self.c3multiplier)
            self.c2slider.setValue(self.dispersionCoefs[1]/self.c2multiplier)
            self.update(True)
        else:
            print 'No global coef file available.'
            
    def setC3(self,intval):
        self.dispersionCoefs[0] = float(intval)*self.c3multiplier
        self.update(True)

    def setC2(self,intval):
        self.dispersionCoefs[1] = float(intval)*self.c2multiplier
        self.update(True)

    def updateInfo(self):
        print self.imageStatsStringBuffer
        outstring = ''
        for item in self.imageStatsStringBuffer:
            if len(item):
                outstring = outstring + item + '\n'
        outstring = outstring + self.getInfoString()
        self.infoLabel.setText(outstring)

    def getInfoString(self,full=True):
        test = self.proc**2.0
        cis = self.processor.computeIntensitySquared(test)
        csnrmax = self.processor.computeSNR(test)
        csnrmed = self.processor.computeSNR(test,mode='median')
        ccm = self.processor.computeContrastMax(test)
        cc = self.processor.computeContrast(test)
        cpc = self.processor.computeProfileContrast(test)
        cmg = self.processor.computeMaxGradients(test)
        chp = self.processor.computeHighPix(test)
        l1 = 'c3=%0.2e \t c2=%0.2e \n'%(self.dispersionCoefs[0],self.dispersionCoefs[1])
        l2 = 'cis=%0.3e\tcsnrmax=%0.3f\tcsnrmed=%0.3f\tccm=%0.3f\tcc=%0.3f\tcpc=%0.3f\tcmg=%0.3e\tchp=%0.3f'%(cis,csnrmax,csnrmed,ccm,cc,cpc,cmg,chp)
        if full:
            return l1+l2
        else:
            return l2
            #return 'c3=%0.2e \t c2=%0.2e \n cis=%0.2e \t'%(,self.LMIN,self.LMAX)

    def toggleLogMode(self):
        self.logMode = not self.logMode
        if self.logMode:
            self.procHigh = 50
            self.procLow = 3
            self.update()
        else:
            self.procHigh = 5000
            self.procLow = 500
            
    def getLogMode(self):
        return self.logMode

    def getRawHigh(self):
        return self.rawHigh/self.limitMultiplier

    def setRawHigh(self,val):
        self.rawHigh = float(val)*self.limitMultiplier
        self.update()

    def getProcHigh(self):
        return self.procHigh/self.limitMultiplier

    def setProcHigh(self,val):
        self.procHigh = float(val)*self.limitMultiplier
        self.update()

    def getRawLow(self):
        return self.rawLow/self.limitMultiplier

    def setRawLow(self,val):
        self.rawLow = float(val)*self.limitMultiplier
        self.update()

    def getProcLow(self):
        return self.procLow/self.limitMultiplier

    def setProcLow(self,val):
        self.procLow = float(val)*self.limitMultiplier
        self.update()


    def checkBox(self,name,setter,getter):
        cb = QtGui.QCheckBox(name)
        initVal = getter()
        if getter():
            cb.toggle()
        cb.stateChanged.connect(setter)
        return cb

    def button(self,name,setter):
        b = QtGui.QPushButton(name)
        b.clicked.connect(setter)
        return b
        

    def spinBox(self,name,setter,getter,increment,aMax=None,aMin=None,rMax=None,rMin=None,suffix=''):
        sb = QtGui.QSpinBox(self)
        sb.setMaximumWidth(100)
        sb.setMaximumHeight(20)
        sbl = QtGui.QLabel(name)
        sbl.setMaximumWidth(100)
        sbl.setMaximumHeight(20)

        initVal = getter()

        maxVal,minVal = None, None
        if aMax is not None and aMin is not None:
            maxVal = aMax
            minVal = aMin
        elif rMax is not None and rMin is not None:
            maxVal = initVal + rMax
            minVal = initVal + rMin
        else:
            maxVal = 2**16
            minVal = -2**16

        sb.setMaximum(maxVal)
        sb.setMinimum(minVal)

        sb.setSingleStep(increment)
        sb.setSuffix(' '+suffix)
        sb.setValue(round(initVal))
        sb.valueChanged[int].connect(setter)

        # do the layout:
        splitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
        splitter.addWidget(sbl)
        splitter.addWidget(sb)
        return splitter


    def regetAndUpdate(self):
        self.update(True)

    def update(self,reget=False):
        print 'update called!!'

        #rcanvasw = self.rcanvas.size().width()
        #rcanvash = self.rcanvas.size().height()
        pcanvasw = self.pcanvas.size().width()
        pcanvash = self.pcanvas.size().height()

        raw = self.getRawImage(reget)
        proc = self.getProcessedImage(reget)

        #self.rqi = QtGui.QImage(raw.data,self.sx,self.sy,self.sx,QtGui.QImage.Format_Indexed8)
        #self.rqi.setColorTable(self._cmap)
        #self.rcanvas.setPixmap(QtGui.QPixmap.fromImage(self.rqi).scaled(rcanvasw,rcanvash,self.aspect))


        self.pqi = QtGui.QImage(proc.data,self.sx,self.sy,self.sx,QtGui.QImage.Format_Indexed8)
        #self.rqi.setColorTable(self._cmap)
        self.pcanvas.setPixmap(QtGui.QPixmap.fromImage(self.pqi).scaled(pcanvasw,pcanvash,self.aspect))
        self.updateInfo()


    def getRawImage(self,reget=False):
        if reget:
            raw = self.processor.getFrame(self.iVol,self.iSlow)
            self.raw = raw

        raw = np.clip(self.raw,self.rawLow,self.rawHigh)
        fval = (raw-raw.min())/(raw.max()-raw.min())
        raw = np.uint8(np.round(fval*255))
        raw = np.clip(raw,0,255)
        raw = raw[::4]
        raw = np.require(raw,np.uint8,'C')
        return raw

    def getProcessedImage(self,reget=False):
        if reget:
            proc = np.abs(self.processor.processFrame(self.raw,dispersionCoefs=self.dispersionCoefs))
            self.proc = proc
            print(self.getInfoString(full=False))
            #print('SNR: %0.2f\t contrast: %0.4f'%(self.processor.computeSNR(proc),self.processor.computeContrast(proc)))
            # self.imageStatsStringBuffer = ['']*5
            self.imageStatsStringBuffer.append(self.getImageStatsString(self.proc))
            self.imageStatsStringBuffer = self.imageStatsStringBuffer[1:]
        if not self.logMode:
            proc = np.clip(self.proc,self.procLow,self.procHigh)
            fval = (proc-proc.min())/(proc.max()-proc.min())
            proc = np.round(fval*255)
            proc = np.clip(proc,0,255)
            proc = np.uint8(proc)
        else:
            proc = np.log10(self.proc)
            proc = np.clip(proc,self.procLow,self.procHigh)
            fval = (proc-proc.min())/(proc.max()-proc.min())
            proc = np.round(fval*255)
            proc = np.clip(proc,0,255)
            proc = np.uint8(proc)
            
        proc = np.require(proc,np.uint8,'C')
        
        return proc
        
    def saveSettings(self):
        np.savetxt(os.path.join('.settings','%s_clims.txt'%self._name),[self.cMin,self.cMax])

    def loadSettings(self):
        vals = np.loadtxt(os.path.join('.settings','%s_clims.txt'%self._name))
        return vals


def main():
    
    app = QtGui.QApplication(sys.argv)

    if len(sys.argv)>1:
        p = OCTTools.Processor(sys.argv[1],True)
        ex = ProcessingWindow(p)
    else:
        ex = ProcessingWindow()

    sys.exit(app.exec_())
    sys.exit()

if __name__ == '__main__':
    main()
