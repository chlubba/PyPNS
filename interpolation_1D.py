"""
Check interpolation performance for simple model scenario


"""

import numpy as np
from scipy.interpolate import Rbf
from scipy.interpolate import NearestNDInterpolator
import matplotlib.pyplot as plt


points = 1000

xObserved = np.squeeze(np.random.rand(1, points))*np.pi*2
yObserved = np.squeeze(np.random.rand(1, points))*np.pi*2

vObserved = np.sin(xObserved)*np.sin(yObserved) # np.add(np.sin(xObserved), np.squeeze(np.random.rand(1, 40))*0.4)

interpolatorN = NearestNDInterpolator(np.transpose(np.vstack([xObserved, yObserved])), vObserved)

interPoints = 100
xInter = yInter = np.linspace(0, np.pi*2, interPoints)
xxInter, yyInter = np.meshgrid(xInter, yInter)

xxInter1D = np.squeeze(xxInter.reshape((1,-1)))
yyInter1D = np.squeeze(yyInter.reshape((1,-1)))

vInter = interpolatorN(np.transpose(np.vstack([xxInter1D, yyInter1D])))

vvInter = vInter.reshape((interPoints,-1))


from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter3D(xObserved, yObserved, vObserved, linewidth=3, color='r')
# ax.plot_surface(xxInter, yyInter, vvInter)
ax.plot_wireframe(xxInter, yyInter, vvInter)
plt.show()


