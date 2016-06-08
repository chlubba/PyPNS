import scipy.ndimage as ndimage
import numpy as np
import matplotlib.pyplot as plt

weights = [1] # ,2,3,4,5,6,5,4,3,2,1]

t = np.tile(np.arange(0,1000,.1), (2, 1))
signalRec = np.sin(2*np.pi/1000*t)

convolved = ndimage.convolve1d(signalRec, weights, axis=-1)

plt.plot(np.transpose(signalRec))
plt.plot(np.transpose(convolved))
plt.show()

