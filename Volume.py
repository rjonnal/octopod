import octopod_config as ocfg



class BScan:
    """A class representing an OCT B-scan.
    """
    
    def __init__(self,fast_start_position,slow_position,depth_start_position,image,labels={}):
        """A class representing an OCT B-scan.
        Args:
            fast_start_position (float): the starting position of the B-scan in the fast dimension, in meters
            slow_position (float): the position of the B-scan in the slow dimension, in meters
            depth_start_position (float): the depth of the top row of pixels, in meters
            image (complex array): the B-scan data
            labels (dict): a dictionary mapping strings (e.g. 'IS/OS') onto one-dimensional arrays of numbers
                such that labels['IS/OS'][10] returns the depth of 'IS/OS' (in pixels) of the A-line with index 10.
        """
        self.image = image
        self.fast_start_position = fast_start_position
        self.slow_position = slow_position
        self.depth_start_position = depth_start_position
        self.labels = labels


        


                 
    
