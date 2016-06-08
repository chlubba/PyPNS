import numpy as np
import time
from PyPN.takeTime import *
from scipy.interpolate import Rbf
from scipy.interpolate import NearestNDInterpolator

points = [7000] # np.power([2],range(10,20))

for point in points:

    print '%i points' % point

    x, y, z, d = np.random.rand(4, point)

    with takeTime('interpolator'):
        coordinates = np.transpose(np.vstack([x, y, z]))
        interpolatorN = NearestNDInterpolator(coordinates, d)
        interpolatorRbf = Rbf(x, y, z, d)  # radial basis function interpolator instance

    xi = yi = zi = np.linspace(0, 1, point)

    with takeTime('interpolation'):
        diN = interpolatorN(xi, yi, zi)   # interpolated values
        diRbf = interpolatorRbf(xi, yi, zi)  # interpolated values

    print np.sqrt(np.mean(np.power((diN - diRbf), 2)))



