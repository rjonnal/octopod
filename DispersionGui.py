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
from time import time,sleep

logging.basicConfig(level=logging.DEBUG)

dz = 10
Z_CUTON = 200+dz
Z_CUTOFF = 350+dz
X_START = 2
N_STEPS = 3

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
        self.coef_3_label.setText('%s / %0.2e'%(self.coef_3_multiplier_text,self.dispersion_coefs[0]))
        self.coef_2_label.setText('%s / %0.2e'%(self.coef_2_multiplier_text,self.dispersion_coefs[1]))

    def change_coefs(self):
        new_c3 = self.coef_3_text.value() * ocfg.dispersion_3_multiplier
        new_c2 = self.coef_2_text.value() * ocfg.dispersion_2_multiplier
        self.dispersion_coefs[0] = new_c3
        self.dispersion_coefs[1] = new_c2

        self.coef_3_max = max(self.coef_3_max,new_c3)
        self.coef_3_min = min(self.coef_3_min,new_c3)
        self.coef_2_max = max(self.coef_2_max,new_c2)
        self.coef_2_min = min(self.coef_2_min,new_c2)
        
        self.show_bscan()


    def zero_coefs(self):
        self.coef_3_text.setValue(0.0)
        self.coef_2_text.setValue(0.0)
        self.change_coefs()
        self.show_bscan()
        
        
    def change_index(self):
        new_index = (self.index_text.value())%(self.raw_vol.shape[0])
        self.index = new_index
        self.show_bscan()

    def optimize(self,obj=np.max,index=0):
        c3_increment = self.coef_3_step_size*self.coef_3_multiplier
        c2_increment = self.coef_2_step_size*self.coef_2_multiplier
        c3_vec = np.arange(self.coef_3_min,self.coef_3_max+c3_increment,c3_increment)
        c2_vec = np.arange(self.coef_2_min,self.coef_2_max+c2_increment,c2_increment)
        scores = np.ones((len(c3_vec),len(c2_vec)))*obj(self.bscan)
        winning_score = -np.inf
        
        c3_winner = self.dispersion_coefs[0]
        c2_winner = self.dispersion_coefs[1]


        for idx3,c3 in enumerate(c3_vec):
            for idx2,c2 in enumerate(c2_vec):
                coefs = [c3,c2,0.0,0.0]
                try:
                    test = self.proc_cache[(index,c3,c2)]
                except KeyError as ke:
                    test = np.abs(process(self.raw_vol[index,:,:],self.k_in,self.k_out,coefs)[Z_CUTON:Z_CUTOFF,X_START:]).T
                    self.proc_cache[(index,c3,c2)] = test
                score = obj(test)
                scores[idx3,idx2] = score
                if score>winning_score:
                    c3_winner = c3
                    c2_winner = c2
                plt.imshow(scores,interpolation='none')
                plt.pause(.1)

        plt.colorbar()
        self.coef_3_text.setValue(c3_winner/self.coef_3_multiplier)
        self.coef_2_text.setValue(c2_winner/self.coef_2_multiplier)
        self.change_coefs()
        self.show_bscan()


    def gradient(self,im):
        return np.mean(np.abs(np.diff(im,axis=0)))
        
    def optimize2(self,obj=None,index=None):
        if obj is None:
            obj = self.gradient
        
        if index is None:
            index = self.index
            
        c3_increment = self.coef_3_step_size*self.coef_3_multiplier
        c2_increment = self.coef_2_step_size*self.coef_2_multiplier

        #c3_vec = np.arange(self.coef_3_min,self.coef_3_max+c3_increment,c3_increment)
        #c2_vec = np.arange(self.coef_2_min,self.coef_2_max+c2_increment,c2_increment)

        c3_vec = np.arange(self.dispersion_coefs[0]-c3_increment*N_STEPS,self.dispersion_coefs[0]+c3_increment*(N_STEPS+1),c3_increment)
        c2_vec = np.arange(self.dispersion_coefs[1]-c2_increment*N_STEPS,self.dispersion_coefs[1]+c2_increment*(N_STEPS+1),c2_increment)
        
        c3_max = np.max(c3_vec)
        c2_max = np.max(c2_vec)
        c3_min = np.min(c3_vec)
        c2_min = np.min(c2_vec)

        c3 = self.dispersion_coefs[0]
        c2 = self.dispersion_coefs[1]

        scores = np.ones((len(c3_vec),len(c2_vec)))*obj(self.bscan)


        print c3_min,c3,c3_max
        print c2_min,c2,c2_max
        while c3_min<=c3<=c3_max and c2_min<=c2<=c2_max:
            print c3_min<=c3<=c3_max
            print c2_min<=c2<=c2_max
            best_score = -np.inf
            for d3 in range(-1,2):
                for d2 in range(-1,2):
                    coefs = [c3+d3*c3_increment,c2+d2*c2_increment,0.0,0.0]
                    try:
                        test = self.proc_cache[(index,c3,c2)]
                    except KeyError as ke:
                        test = np.abs(process(self.raw_vol[index,:,:],self.k_in,self.k_out,coefs)[Z_CUTON:Z_CUTOFF,X_START:]).T
                        self.proc_cache[(index,c3,c2)] = test
                    score = obj(test)
                    if score>best_score:
                        best_d3 = d3
                        best_d2 = d2
                        best_score = score
            c3 = c3+best_d3*c3_increment
            c2 = c2+best_d2*c2_increment
            i3 = np.argmin(np.abs(c3-c3_vec))
            i2 = np.argmin(np.abs(c2-c2_vec))
            scores[i3,i2] = best_score

            self.canvas.setImage(scores)
            
            if best_d3==0 and best_d2==0:
                break

        c3_winner = c3
        c2_winner = c2
        
        self.coef_3_text.setValue(c3_winner/self.coef_3_multiplier)
        self.coef_2_text.setValue(c2_winner/self.coef_2_multiplier)
        self.change_coefs()
        self.show_bscan()

    def compute_stats(self):
        try:
            self.bscan_max = np.max(self.bscan)
            self.bscan_min = np.min(self.bscan)
            self.bscan_mean = np.mean(self.bscan)
            self.bscan_med = np.median(self.bscan)
            self.bscan_var = np.var(self.bscan)
            self.bscan_range = self.bscan_max - self.bscan_min
        except Exception as e:
            self.bscan_max = 0.0
            self.bscan_min = 0.0  
            self.bscan_mean = 0.0 
            self.bscan_med = 0.0  
            self.bscan_var = 0.0  
            self.bscan_range = 0.0
        self.stat_text = 'max = %0.3e\nmin = %0.3e\nmean=%0.3e\nmed=%0.3e\nvar=%0.3e\nrange=%0.3e'%(self.bscan_max,self.bscan_min,self.bscan_mean,self.bscan_med,self.bscan_var,self.bscan_range)
        self.stats_label.setText(self.stat_text)
            
    def init_UI(self):

        self.coef_3_multiplier = ocfg.dispersion_3_multiplier
        self.coef_2_multiplier = ocfg.dispersion_2_multiplier
        self.coef_3_step_size = ocfg.dispersion_3_step_size
        self.coef_2_step_size = ocfg.dispersion_2_step_size
        self.coef_3_multiplier_text = 'x 10^%d'%(np.log10(self.coef_3_multiplier))
        self.coef_2_multiplier_text = 'x 10^%d'%(np.log10(self.coef_2_multiplier))

        ## Create some widgets to be placed inside
        btn_open = QtGui.QPushButton('&open file')
        btn_quit = QtGui.QPushButton('&quit')
        btn_write = QtGui.QPushButton('&write to hdf5')
        btn_optimize = QtGui.QPushButton('o&ptimize')
        btn_clear_cache = QtGui.QPushButton('&clear cache')
        btn_zero_coefs = QtGui.QPushButton('&zero coefs')
        
        btn_open.clicked.connect(self.open_file)
        btn_quit.clicked.connect(self.close)
        btn_write.clicked.connect(self.write_coefs)
        btn_optimize.clicked.connect(self.optimize2)
        btn_clear_cache.clicked.connect(self.clear_cache)
        btn_zero_coefs.clicked.connect(self.zero_coefs)
        
        self.coef_3_text = QtGui.QDoubleSpinBox()
        self.coef_2_text = QtGui.QDoubleSpinBox()

        self.index_text = QtGui.QSpinBox()
        self.index_text.setMaximum(1000)
        self.index_text.setMinimum(0)

        

        self.coef_3_text.setMaximum(ocfg.dispersion_3_max)
        self.coef_3_text.setMinimum(ocfg.dispersion_3_min)
        self.coef_3_text.setSingleStep(self.coef_3_step_size)
        
        self.coef_2_text.setMaximum(ocfg.dispersion_2_max)
        self.coef_2_text.setMinimum(ocfg.dispersion_2_min)
        self.coef_2_text.setSingleStep(self.coef_2_step_size)

        self.coef_3_label = QtGui.QLabel('%0.2e'%self.dispersion_coefs[0])
        self.coef_2_label = QtGui.QLabel('%0.2e'%self.dispersion_coefs[1])


        self.stats_label = QtGui.QLabel('')
        
        self.canvas = pg.ImageView()

        self.statusbar = QtGui.QStatusBar()
        
        ## Create a grid layout to manage the widgets size and position
        layout = QtGui.QGridLayout()
        layout.setSpacing(10)
        self.setLayout(layout)

        ## Add widgets to the layout in their proper positions
        layout.addWidget(self.canvas, 0, 0, 5, 9)  # plot goes on right side, spanning 3 rows
        layout.addWidget(self.stats_label, 0, 9, 5, 1)
        layout.addWidget(btn_open, 5, 0, 1, 1)   # button goes in upper-left
        layout.addWidget(btn_optimize, 6, 0, 1, 1)
        layout.addWidget(btn_write, 7, 0, 1, 1)   # text edit goes in middle-left
        layout.addWidget(btn_quit, 8, 0, 1, 1)   # text edit goes in middle-left
        layout.addWidget(btn_clear_cache, 5, 9, 1, 1)
        layout.addWidget(btn_zero_coefs, 6, 9, 1, 1)
        layout.addWidget(self.statusbar,9,0,1,10)
        
        label_3 = QtGui.QLabel('&3rd Order')
        label_2 = QtGui.QLabel('&2nd Order') 
        label_i = QtGui.QLabel('Frame &index')        
        label_3.setBuddy(self.coef_3_text)
        label_2.setBuddy(self.coef_2_text)
        label_i.setBuddy(self.index_text)
        
        layout.addWidget(label_3, 5, 1, 1, 1)
        layout.addWidget(label_2, 6, 1, 1, 1)
        layout.addWidget(label_i, 7, 1, 1, 1)
        
        layout.addWidget(self.coef_3_text, 5, 2, 1, 1)  # list widget goes in bottom-left
        layout.addWidget(self.coef_2_text, 6, 2, 1, 1)  # list widget goes in bottom-left
        layout.addWidget(self.index_text, 7, 2, 1, 1)

        layout.addWidget(self.coef_3_label, 5, 3, 1, 1)
        layout.addWidget(self.coef_2_label, 6, 3, 1, 1)

        self.setMinimumWidth(1800)
        self.setMinimumHeight(1000)

        self.compute_stats()
        
        self.show()


    def open_file(self,fname=None):
        
        if fname is None:
            try:
                fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file',self.working_dir)[0]
            except:
                fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file',ocfg.data_root)[0]
                
        if os.path.exists(fname):
            self.working_dir,junk = os.path.split(fname)
            self.index = 10
            self.statusbar.showMessage('Working file: %s'%fname)
            self.h5 = H5(fname)
            self.dispersion_optimizer = DispersionOptimizer(fname)
            self.raw_vol = self.h5.get('raw_data')[0,:,:,:]
            self.k_in = self.h5.get('k_in')
            self.k_out = self.h5.get('k_out')
            self.dispersion_coefs = self.h5.get('dispersion/coefficients')
            
            self.coef_3_min = self.dispersion_coefs[0]
            self.coef_3_max = self.dispersion_coefs[0]
            self.coef_2_min = self.dispersion_coefs[1]
            self.coef_2_max = self.dispersion_coefs[1]

            self.coef_3_text.setValue(self.dispersion_coefs[0]/self.coef_3_multiplier)
            self.coef_2_text.setValue(self.dispersion_coefs[1]/self.coef_2_multiplier)
            self.index_text.setValue(self.index)
            

            self.coef_3_text.valueChanged.connect(self.change_coefs)
            self.coef_2_text.valueChanged.connect(self.change_coefs)
            self.index_text.valueChanged.connect(self.change_index)
            
            self.proc_cache = {}
            self.show_bscan()
            
    def write_coefs(self):
        old_coefs = self.h5.get('dispersion/coefficients')
        self.h5.put('dispersion/coefficients',self.dispersion_coefs)
        self.h5.put('dispersion/old_coefficients',old_coefs)

        ddb = H5(ocfg.dispersion_database)
        did = self.h5.get('IDs/dataset_id')[()]
        did_key = '%d'%did
        ddb.put(did_key,self.dispersion_coefs)
        

    def show_bscan(self,index=None):
        if index is None:
            index = self.index
        try:
            self.bscan = self.proc_cache[(index,self.dispersion_coefs[0],self.dispersion_coefs[1])]
        except KeyError as ke:
            self.bscan = np.abs(process(self.raw_vol[index,:,:],self.k_in,self.k_out,self.dispersion_coefs)[Z_CUTON:Z_CUTOFF,X_START:]).T
            self.proc_cache[(index,self.dispersion_coefs[0],self.dispersion_coefs[1])] = self.bscan
        self.compute_stats()
        self.set_coef_labels()
        self.canvas.setImage(self.bscan)
        #print sys.getsizeof(self.proc_cache)

    def clear_cache(self):
        self.proc_cache = {}


def main():
    ## Always start by initializing Qt (only once per application)
    app = QtGui.QApplication([])
    ex = Window()

    try:
        ex.open_file(sys.argv[1])
    except:
        print 'No file specified.'
    
    ## Start the Qt event loop
    app.exec_()


if __name__=='__main__':
    main()
