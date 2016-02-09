from PySide import QtGui  # (the example applies equally well to PySide)
import pyqtgraph as pg
import sys,os
from octopod import *
import numpy as np
import octopod_config as ocfg
from octopod.Processor import process
from octopod import DispersionOptimizer
from pyqtgraph import ImageView




class Window(QtGui.QWidget):
    def __init__(self):
        super(Window,self).__init__()
        self.init_UI()

    def init_UI(self):


## Define a top-level widget to hold everything
#w = QtGui.QWidget()

        ## Create some widgets to be placed inside
        btn_open = QtGui.QPushButton('open file')
        btn_quit = QtGui.QPushButton('quit')
        coef_3_text = QtGui.QDoubleSpinBox()
        coef_2_text = QtGui.QDoubleSpinBox()

        coef_3_text.setMaximum(ocfg.DISPERSION_3_MAX)
        coef_3_text.setMinimum(ocfg.DISPERSION_3_MIN)
        coef_2_text.setMaximum(ocfg.DISPERSION_2_MAX)
        coef_2_text.setMinimum(ocfg.DISPERSION_2_MIN)

        
        plot = pg.ImageView()

        ## Create a grid layout to manage the widgets size and position
        layout = QtGui.QGridLayout()
        layout.setSpacing(10)
        self.setLayout(layout)

        ## Add widgets to the layout in their proper positions
        layout.addWidget(plot, 0, 0, 5, 10)  # plot goes on right side, spanning 3 rows
        layout.addWidget(btn_open, 5, 0, 1, 1)   # button goes in upper-left
        layout.addWidget(btn_quit, 6, 0, 1, 1)   # text edit goes in middle-left
        layout.addWidget(coef_3_text, 5, 1, 1, 1)  # list widget goes in bottom-left


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
        return process(self.raw_vol[index,:,:],self.k_in,self.k_out,self.coefficients)[:,2:]



def main():
    ## Always start by initializing Qt (only once per application)
    app = QtGui.QApplication([])

    ex = Window()

    ## Start the Qt event loop
    app.exec_()


if __name__=='__main__':
    main()
