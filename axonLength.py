import numpy as np
import PyPN
import matplotlib.pyplot as plt

lengthOfBundle = 1000

# bundle guide
segmentLengthAxon = 30
bundleGuide = PyPN.createGeometry.get_bundle_guide_straight_radius(lengthOfBundle, segmentLengthAxon, radius=200)

RDCs = np.arange(0, .9, 0.1)
lengthAxon = []
for RDC in RDCs:
    axonCoords = PyPN.createGeometry.create_random_axon(bundleGuide, [0,0], segmentLengthAxon,
                                               randomDirectionComponent=RDC)

    lengthAxon.append(PyPN.createGeometry.length_from_coords(axonCoords))

plt.plot(RDCs, np.array(lengthAxon)/lengthAxon[0])
plt.show()



