import numpy as np

nx, ny = (3, 2)
x = np.linspace(0, 1, nx)
y = np.linspace(0, 1, ny)
xv, yv = np.meshgrid(x, y)

xvLin = xv.reshape([1,-1])
yvLin = yv.reshape([1,-1])

print xvLin
print yvLin