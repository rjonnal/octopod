import matplotlib
#matplotlib.use('Qt5Agg')
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
import numpy as np
import sys,os
import h5py

from ActivityDetector import *
from AcquisitionParameterFile import *
from Dataset import *
from DatabaseManager import *
from Misc import *
from DataStore import *
from FileManager import *
from Processor import *
from DispersionOptimizer import *
from Model import *
import utils
#from movie import *
from VisualizationTools import *
from Reporter import *
from Series import *
from Cropper import *
from Flipper import *
from StackAverager import *
#from StripRegistrar import *
from Cones import *
from BScanMaker import *
from DispersionUnifier import *
