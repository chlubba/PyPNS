import PyPN
import numpy as np
import matplotlib.pyplot as plt

# set length of bundle and number of axons
lengthOfBundle = 400000
segmentLengthAxon = 30
bundleGuide = PyPN.createGeometry.get_bundle_guide_straight(lengthOfBundle, segmentLengthAxon)

extracellulars = []
extracellulars.append(PyPN.Extracellular.homogeneous(sigma=1.2))
extracellulars.append(PyPN.Extracellular.precomputedFEM(bundleGuide))
extracellulars.append(PyPN.Extracellular.interpolator(bundleGuide))

zPositions = np.arange(-15000, 15000, 10)
sourceCurrents = np.ones(1)

colors = np.array(((0.,0.,0.), (230., 159., 0.), (86., 180., 233.), (0., 158., 115.)))/255
lineStyles = ['-', '--', '-.']

legends = ['homogeneous', 'radial inhomogeneous', 'cuff']
for extraInd, extracellular in enumerate(extracellulars):

    receiverPositions = np.hstack((zPositions[:,np.newaxis], np.ones((len(zPositions),1))*235, np.zeros((len(zPositions),1))))

    for xPInd, xP in enumerate([0, 180]): # [0, 90, 180]

        if extraInd == 2:
            xP = -xP
        sourcePositions = np.array([0, 0, xP])

        v = extracellular.calculate_extracellular_potential(np.transpose(sourcePositions[:, np.newaxis]), sourceCurrents[:, np.newaxis], receiverPositions)
        plt.plot(zPositions, v/np.max(v), lineStyles[xPInd], label=legends[extraInd], color=colors[extraInd])

plt.legend()

import os
xlimits = ([-15000,15000], [-500,500])
figureNames = ['profileZFull.eps', 'profuleZzoomed.eps']
for limInd in range(2):
    plt.xlim(xlimits[limInd])

    plt.savefig(os.path.join('/home/carl/Dropbox/_Exchange/Project/PyPN Paper/PythonFigureOutput', figureNames[limInd]),
            format='eps', dpi=300)

plt.show()

# sourcePositions = np.hstack((np.zeros((len(zPositions),2)), zPositions[:,np.newaxis]))
# receiverPositions = [235, 0, 0]
# homo = PyPN.Extracellular.homogeneous
# v = homo.calculate_extracellular_potential(sourcePositions, sourceCurrents, receiverPositions)
