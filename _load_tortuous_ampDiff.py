import numpy as np
import matplotlib.pyplot as plt

amplitudeMatrix = np.load('/media/carl/4ECC-1C44/PyPN/tortuous/amplitudeDiff.npy')

print amplitudeMatrix

RDCs = [0, 0.2, 0.4, 0.6, 0.8, 1.]
elecRadii = [50, 100, 150, 200, 250, 300]

RDCsMat, elecRadiiMat = np.meshgrid(RDCs, elecRadii)

CS = plt.contour(elecRadii, RDCs, amplitudeMatrix)
plt.clabel(CS, inline=1, fontsize=10)
plt.xlabel('electrode radius [$\mu$m]')
plt.ylabel('random direction component $\\alpha$')
plt.show()

