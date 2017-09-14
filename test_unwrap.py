from octopod import utils
import numpy as np
from matplotlib import pyplot as plt

phases = np.random.randn(100)*.2
phases = phases%(2*np.pi)
plt.subplot(1,2,1)
plt.plot(phases,'ks')
print np.var(phases)
uphases = utils.unwrap(phases)
plt.subplot(1,2,2)
plt.plot(uphases,'ks')
print np.var(uphases)
plt.show()
