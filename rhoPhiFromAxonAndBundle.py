import numpy as np
from math import hypot
# import time
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# problem: need to know the bundle-guide segment associated to the axon

# t0 = time.time()

bundleStartPoint = np.array([0, 0, 0])
bundleDirection = np.array([1, 0, 0])
bundleSementLength = 200 # um

axonSegmentPosition = np.array([100, 50, 50])

electrodePosition = np.array([0,1000, 0])

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# points = np.vstack([bundleStartPoint, axonSegmentPosition])
# ax.scatter(points[:,0], points[:,1], points[:,2])
# plt.show()

# for i in range(100000):
#
#     # compontent parallel to bundle direction
#     elecPosPar = np.inner(electrodePosition - bundleStartPoint, bundleDirection)
#     axonPosPar = np.inner(axonSegmentPosition - bundleStartPoint, bundleDirection)
#
#     # normal from bundle guide to axon position
#     axonDir2D = axonSegmentPosition - (bundleStartPoint + axonPosPar*bundleDirection)
#     axonDir2DNorm = axonDir2D/np.linalg.norm(axonDir2D)
#
#     # # 2nd vector of new coordinate system
#     # axonDir2DNorm2 = np.cross(axonDir2DNorm, bundleDirection)
#
#     # electrode shifted to origin
#     electrodeVector = electrodePosition - (bundleStartPoint + elecPosPar * bundleDirection)
#
#     # electrode coordinates projected onto new base vectors
#     elecX = np.inner(electrodeVector, axonDir2DNorm)
#
#     elecYVec = electrodeVector - elecX * axonDir2DNorm
#     elecY = np.linalg.norm(elecYVec)
#
#     # new coordinates
#     r = hypot(elecX, elecY)
#     phi = np.arctan2(elecY, elecX)
#     z = elecPosPar - axonPosPar

for i in range(100000):

    # ----- axon segment part -----------

    # compontent parallel to bundle direction
    axonPosPar = np.inner(axonSegmentPosition - bundleStartPoint, bundleDirection)

    # normal from bundle guide to axon position
    axonDir2D = axonSegmentPosition - (bundleStartPoint + axonPosPar*bundleDirection)
    axonDir2DNorm = axonDir2D/np.linalg.norm(axonDir2D)

    # ----- electrode position part -----

    # compontent parallel to bundle direction
    elecPosPar = np.inner(electrodePosition - bundleStartPoint, bundleDirection)

    # electrode shifted to origin
    electrodeVector = electrodePosition - (bundleStartPoint + elecPosPar * bundleDirection)

    # electrode coordinates projected onto new base vectors
    elecX = np.inner(electrodeVector, axonDir2DNorm)

    elecYVec = electrodeVector - elecX * axonDir2DNorm
    elecY = np.linalg.norm(elecYVec)

    # new coordinates
    r = hypot(elecX, elecY)
    phi = np.arctan2(elecY, elecX)
    z = elecPosPar - axonPosPar

# for bs in bundleSegments:
#     for a in axonsOfBundleSegment:
#         calculate axonDir2DNorms
#     for e in electrodes calculate:
#         for a in axonsOfBundleSegment:
#               calculate electrode position and interpolate potential

# print time.time() - t0

print elecX
print elecY

# phi = np.arctan2()
