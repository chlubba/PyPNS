import numpy as np
import os

timeSteps = 10
path = '/media/carl/4ECC-1C44/PyPN/voltageFieldDummy2'

for timeStep in range(timeSteps):

    x_ = y_ = z_ = np.linspace(0, 2000, 10)

    xx, yy, zz = np.meshgrid(x_, y_, z_)

    x = np.squeeze(xx.reshape(1,-1))
    y = np.squeeze(yy.reshape(1, -1))
    z = np.squeeze(zz.reshape(1, -1))

    v = np.sin(2*np.pi/2000*x + 2*np.pi/300*timeStep) * np.sin(2*np.pi/2000*y + 2*np.pi/300*timeStep) * \
        np.sin(2 * np.pi / 2000 * z + 2 * np.pi / 300 * timeStep)

    data = np.vstack([x,y,z,v])

    filename = 'voltageField'+str(timeStep)+'.dat'
    location = os.path.join(path, filename)

    np.savetxt(location, data)