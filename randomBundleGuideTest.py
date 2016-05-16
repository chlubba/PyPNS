import PyPN
import numpy as np
import matplotlib.pyplot as plt

initCoords = np.array([[0, 0, 0], [50, 0, 0], [10000, 0, 0]])

# bundleGuide = PyPN.createGeometry.create_random_axon(initCoords, 10000, [0,0], 100, randomDirectionComponent= 0.5 )

segmentLength = 200

numPoints = 20
# randomDelta = np.random.uniform(-1,1,(numPoints, 3))

randomDeltaX = np.random.uniform(0,1,numPoints)
randomDeltaYZ = np.random.uniform(-1,1,(numPoints,2))

randomDelta = np.column_stack((randomDeltaX, randomDeltaYZ))

for i in range(numPoints):
    randomDelta[i,:] = randomDelta[i,:]/np.linalg.norm(randomDelta[i,:])

bundleGuide = np.cumsum(randomDelta,0)

bundleGuideScaled = bundleGuide*segmentLength

fig = plt.figure()
ax = fig.gca(projection='3d')

X = bundleGuideScaled[:,0]
Y = bundleGuideScaled[:,1]
Z = bundleGuideScaled[:,2]

ax.plot(X, Y, Y) # label='axon '+str(axonID),

max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max() / 2.0

mid_x = (X.max()+X.min()) * 0.5
mid_y = (Y.max()+Y.min()) * 0.5
mid_z = (Z.max()+Z.min()) * 0.5
ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)

plt.show()