import cPickle as pickle
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm
import matplotlib.colors as colors
import numpy as np

condVelDict = pickle.load(open(os.path.join('/media/carl/4ECC-1C44/PyPN/condVel', 'conductionVelocitiesMyelinated.dict'), "rb" ))

axonType = condVelDict['axonType']
diameters = condVelDict['diameters']
temperatures = condVelDict['temperatures']
velocityArray = condVelDict['velocityArray']



jet = plt.get_cmap('jet')
cNorm = colors.Normalize(vmin=0, vmax=len(temperatures) - 1)
scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

for temperatureInd, temperature in enumerate(temperatures):
    vAPs = velocityArray[temperatureInd]


    colorVal = scalarMap.to_rgba(temperatureInd)
    plt.plot(diameters, vAPs, label=str(temperature) + ' $^\circ$C', color=colorVal)

# # save the bundle to disk
# PyPN.save_bundle(bundle)

plt.plot([0, 4.7], [0, 23.5], linestyle='--', color=np.array([1, 1, 1])*0.7, label='theory')

plt.xlabel('diameter [um]')
plt.ylabel('conduction velocity [m/s]')
plt.title('Myelinated Axon')
plt.legend(loc='best')
plt.grid()

plt.show()

