import octopod_config as ocfg

class ALine:
    """A class representing single OCT A-lines.
    """
    def __init__(self,a,depth_sampling_interval,image_fast_coordinate,image_slow_coordinate,labels={}):
        """A class representing single OCT A-lines.

        Args:
            a (a list or array of floats): this represents the complex OCT signal at the corresponding depths z.
        """
        self.z = z
        self.amplitude = np.abs(a)
        self.phase = np.angle(a)
        self.labels = labels


class Volume:

    
