import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt

dz = 0.000001
z = np.arange(-0.02, 0.02, dz)
a = 10.
b = 0.01

a_c = a/100
b_c = b/4.7
exp_c = 11

smoothWidth = 0.002
smoothSamples = smoothWidth/dz

def fn (z):
    return np.maximum(0, a*(np.subtract(1,np.abs(np.divide(z,b)))))#  + np.maximum(0, a/10*(np.subtract(1,np.abs(np.divide((z-b),b))))**3) #

def fn_corner_r (z):
    return np.maximum(0, a_c*(np.subtract(1,np.abs(np.divide((z-b),b_c))))**exp_c) # np.maximum(0, a*(np.subtract(1,np.abs(np.divide(z,b))))) +

def fn_corner_l (z):
    return np.maximum(0, a_c*(np.subtract(1,np.abs(np.divide((z+b),b_c))))**exp_c) # np.maximum(0, a*(np.subtract(1,np.abs(np.divide(z,b))))) +

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


sharpFn = fn(z)
smoothedFn = smooth(fn(z),smoothSamples)
combinedFn = np.concatenate([smoothedFn[z < -b + smoothWidth], sharpFn[np.logical_and(z >= -b + smoothWidth, z <= b - smoothWidth)], smoothedFn[z > b - smoothWidth]])

sharpFnOneSide = np.maximum(0, a*(np.add(1,np.divide(z,b))))
smoothedFnOneSide = smooth(sharpFnOneSide, smoothSamples)
smoothedFnOneSideToMiddle = smoothedFnOneSide[0:np.shape(smoothedFnOneSide)[0]/2]
smoothedFnTwoSides = np.concatenate([smoothedFnOneSideToMiddle, np.fliplr([smoothedFnOneSideToMiddle])[0]])
# plt.plot(sharpFnOneSide)

# plt.plot(smoothedFnTwoSides)
# plt.plot(fn(z))
# plt.plot(fn_corner_r(z))
# plt.plot(fn_corner_l(z))
# plt.plot(np.diff(np.diff(fn(z) + fn_corner_r(z) + fn_corner_l(z))))
# plt.plot(smooth(fn(z),1000))
plt.plot(np.diff(np.diff(smoothedFnTwoSides)))
# plt.plot(smooth(fn(z),10) - fn(z))
plt.show()