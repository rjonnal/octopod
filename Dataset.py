import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import logging

class Dataset:

    def __init__(self,raw_data_filename):
        """Initiate a Dataset object based on a .unp raw OCT data file."""
        self.raw_data_filename = raw_data_filename
        
        self.logger = logging.getLogger(__name__)
        self.logger.info('Creating Dataset based on %s.'%raw_data_filename)

        
        
