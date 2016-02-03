import sys,os
from PySide import QtGui, QtCore
from octopod import *
import numpy as np
from matplotlib import pyplot as plt
import octopod_config as ocfg
from octopod.Processor import process
from octopod import DispersionOptimizer

class DispersionGui(QtGui.QMainWindow):

    def __init__(self):
        super(DispersionGui, self).__init__()
        self.init_UI()
        
    def init_UI(self):

        # set up menu actions

        open_action = QtGui.QAction(QtGui.QIcon('open.png'), '&Open', self)
        open_action.setShortcut('Ctrl+O')
        open_action.setStatusTip('Open a file')
        open_action.triggered.connect(self.open_file)
        
        exit_action = QtGui.QAction(QtGui.QIcon('exit.png'), '&Exit', self)
        exit_action.setShortcut('Ctrl+Q')
        exit_action.setStatusTip('Exit application')
        exit_action.triggered.connect(self.close)

        self.statusBar()
        
        menubar = self.menuBar()
        file_menu = menubar.addMenu('&File')
        file_menu.addAction(open_action)
        file_menu.addAction(exit_action)
        
        #self.setGeometry(300, 300, 250, 150)
        #self.setWindowTitle('Menubar')
        self.setMinimumWidth(1800)
        self.setMinimumHeight(1000)
        self.show()

    def open_file(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file',ocfg.data_root)[0]
        if os.path.exists(fname):
            self.setStatusTip('Working file: %s'%fname)
            self.h5 = h5py.File(fname)
            self.dispersion_optimizer = DispersionOptimizer(self.h5)
            


def main():
    
    app = QtGui.QApplication(sys.argv)

    ex = DispersionGui()

    sys.exit(app.exec_())
    sys.exit()

if __name__ == '__main__':
    main()
