import cPickle as pickle
import os
from matplotlib import pyplot as plt
import numpy as np

saveLocation = '/media/carl/4ECC-1C44/PyPN'

interFileName = 'conductionVelocitiesInter.dict'
origFileName  = 'conductionVelocitiesOrig.dict'

conVelInter = pickle.load(open(os.path.join(saveLocation, interFileName), "rb" ))
conVelOrig = pickle.load(open(os.path.join(saveLocation, origFileName), "rb" ))

colors = ('blue', 'red')
labels = ('original', 'interpolated')

counter = 0
for conVel in (conVelOrig, conVelInter):

    myel = conVel['myel']

    diameters = myel['diams']

    if counter == 0:
        # find diameter closest to give nones in McIntyre's paper
        possibleDiams = np.array([5.7, 7.3, 8.7, 10., 11.5, 12.8, 14., 15., 16.])

        diams = []
        for diameter in diameters:

            diffArray = np.abs(np.subtract(possibleDiams, diameter))
            diameterIndex = diffArray.argmin()

            fiberD = possibleDiams[diameterIndex]

            diams.append(fiberD)
        diams = np.array(diams)
    else:
        diams = diameters

    velocityMean = myel['meanVelocity']
    velocityStd = myel['stdVelocity']

    plt.scatter(diams, velocityMean, color=colors[counter], label=labels[counter])
    plt.errorbar(diams, velocityMean, yerr=velocityStd, linestyle='None')

    counter += 1

plt.legend()
plt.show()