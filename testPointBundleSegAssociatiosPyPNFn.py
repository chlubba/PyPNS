import numpy as np
import PyPN
import PyPN.createGeometry as createGeometry
from PyPN.takeTime import takeTime
import PyPN.extracellularBackend as eB

# set length of bundle and number of axons
lengthOfBundle = 40000

# bundle guide
segmentLengthAxon = 30
# bundleGuide = np.array([[0,0,0,200], [200,0,0,200], [400,100,0,200], [600,0,0,200], [800,0,0,200]])
bundleGuide = np.array(PyPN.createGeometry.get_bundle_guide_straight_radius(lengthOfBundle, segmentLengthAxon))
# bundleGuide = np.array(PyPN.createGeometry.get_bundle_guide_random_radius(lengthOfBundle, segmentLengthAxon))

# calculate the random axon coordinates
axonCoords = np.array(createGeometry.create_random_axon(bundleGuide, [0,1], segmentLengthAxon, randomDirectionComponent=0.6))
numAxonSegments = np.shape(axonCoords)[0]

bundleGuideLength = np.shape(bundleGuide)[0]
# print bundleGuideLength

bundleGuideStripped = bundleGuide[:,0:3] # no radius



# with takeTime('calculating squared distances'):
#     segmentAssociations = eB.associatePointToBundleSegs(axonCoords, bundleGuide)
#     distPerpendiculars, lengthAlongs, angles = eB.spatialRelation(axonCoords, bundleGuide, segmentAssociations)
#
# print 'mean %3.3f std %3.3f' % (np.mean(distPerpendiculars), np.std(distPerpendiculars))


electrodePos = np.array([[100, 0, 100], [100, 0, -100]])

segmentAssociationsRec = eB.associatePointToBundleSegs(electrodePos, bundleGuide)
distPerpendicularRec, lengthAlongRec, anglesRec = eB.spatialRelation(electrodePos, bundleGuide, segmentAssociationsRec)

segmentAssociationsSource = eB.associatePointToBundleSegs(axonCoords, bundleGuide)
distPerpendicularsSource, lengthAlongsSource, anglesSource = eB.spatialRelation(axonCoords, bundleGuide, segmentAssociationsSource)

# loop over all recording electrodes
for recInd in range(np.shape(distPerpendicularRec)[0]):

    distPerpendicularRecTemp = distPerpendicularRec[recInd]
    lengthAlongRecTemp = lengthAlongRec[recInd]
    angleRecTemp = anglesRec[recInd]

    # distance between source and recording positions
    distAlongTemp = lengthAlongsSource - lengthAlongRecTemp

    # angle between them
    anglesTemp = anglesSource - angleRecTemp