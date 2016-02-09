from PySide import QtGui  # (the example applies equally well to PySide)
import pyqtgraph as pg
import sys,os
from octopod import *
import numpy as np
import octopod_config as ocfg
from octopod.Processor import process
from octopod import DispersionOptimizer
from pyqtgraph import ImageView
import logging

logging.basicConfig(level=logging.DEBUG)

class Window(QtGui.QWidget):
    def __init__(self):
        super(Window,self).__init__()
        self.dispersion_coefs = [0.0,0.0,0.0,0.0]
        self.init_UI()
        self.logger = logging.getLogger(__name__)

    def report(self):
        poly = '%0.1e x^3 + %0.1e x^2'%(self.coef_3_text.value(),self.coef_2_text.value())
        self.logger.info(poly)

    def set_coef_labels(self):
        self.coef_3_label.setText('%0.2e'%(self.dispersion_coefs[0]))
        self.coef_2_label.setText('%0.2e'%(self.dispersion_coefs[1]))

    def change_coefs(self,newval):
        self.dispersion_coefs[0] = self.coef_3_text.value() * ocfg.dispersion_3_multiplier
        self.dispersion_coefs[1] = self.coef_2_text.value() * ocfg.dispersion_2_multiplier
        self.bscan = self.get_bscan()
        
    def init_UI(self):

        ## Create some widgets to be placed inside
        btn_open = QtGui.QPushButton('open file')
        btn_quit = QtGui.QPushButton('quit')

        btn_open.clicked.connect(self.open_file)
        btn_quit.clicked.connect(self.report)
        
        self.coef_3_text = QtGui.QDoubleSpinBox()
        self.coef_2_text = QtGui.QDoubleSpinBox()

        self.coef_3_multiplier = ocfg.dispersion_3_multiplier
        self.coef_2_multiplier = ocfg.dispersion_2_multiplier

        self.coef_3_text.setMaximum(ocfg.dispersion_3_max)
        self.coef_3_text.setMinimum(ocfg.dispersion_3_min)
        self.coef_3_text.setSingleStep(.01)
        
        self.coef_2_text.setMaximum(ocfg.dispersion_2_max)
        self.coef_2_text.setMinimum(ocfg.dispersion_2_min)
        self.coef_2_text.setSingleStep(.01)

        self.coef_3_label = QtGui.QLabel('%0.2e'%self.dispersion_coefs[0])
        self.coef_2_label = QtGui.QLabel('%0.2e'%self.dispersion_coefs[1])

        self.coef_3_text.valueChanged.connect(self.change_coefs)
        self.coef_2_text.valueChanged.connect(self.change_coefs)
        
        
        plot = pg.ImageView()

        ## Create a grid layout to manage the widgets size and position
        layout = QtGui.QGridLayout()
        layout.setSpacing(10)
        self.setLayout(layout)

        ## Add widgets to the layout in their proper positions
        layout.addWidget(plot, 0, 0, 5, 10)  # plot goes on right side, spanning 3 rows
        layout.addWidget(btn_open, 5, 0, 1, 1)   # button goes in upper-left
        layout.addWidget(btn_quit, 6, 0, 1, 1)   # text edit goes in middle-left
        layout.addWidget(self.coef_3_text, 5, 1, 1, 1)  # list widget goes in bottom-left
        layout.addWidget(self.coef_2_text, 6, 1, 1, 1)  # list widget goes in bottom-left
        layout.addWidget(self.coef_3_label, 5, 2, 1, 1)
        layout.addWidget(self.coef_2_label, 6, 2, 1, 1)

        self.setMinimumWidth(1800)
        self.setMinimumHeight(1000)
        
        self.show()


    def open_file(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file',ocfg.data_root)[0]
        if os.path.exists(fname):
            self.setStatusTip('Working file: %s'%fname)
            self.h5 = h5py.File(fname)
            self.dispersion_optimizer = DispersionOptimizer(self.h5)
            self.raw_vol = self.h5['raw_data'][0,:,:,:]
            self.k_in = self.h5['k_in'][:]
            self.k_out = self.h5['k_out'][:]
            self.coefficients = self.h5['dispersion']['coefficients'][:]
            self.bscan = self.get_bscan()

    def get_bscan(self,index=0):
        return process(self.raw_vol[index,:,:],self.k_in,self.k_out,self.dispersion_coefs)[:,2:]



def main():
    ## Always start by initializing Qt (only once per application)
    app = QtGui.QApplication([])

    ex = Window()

    ## Start the Qt event loop
    app.exec_()


if __name__=='__main__':
    main()
