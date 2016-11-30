import numpy as np
import matplotlib.pyplot as plt

# maxcorMatrix1 = np.load('/media/carl/4ECC-1C44/PyPN/tortuous/tortuous.dat.npy')
# maxcorMatrix2 = np.load('/media/carl/4ECC-1C44/PyPN/tortuous/tortuous2.npy')

# maxcorMatrix1 = np.load('/media/carl/4ECC-1C44/PyPN/tortuous/tortuous_50_300_maxcorr_box100.npy')
# maxcorMatrix2 = np.load('/media/carl/4ECC-1C44/PyPN/tortuous/tortuous_400_1000_maxcorr_box100.npy')

maxcorMatrix1 = np.load('/media/carl/4ECC-1C44/PyPN/tortuous/tortuous_50_300_maxcorr_box100_10_runs.npy')

# maxcorMatrix = np.column_stack([maxcorMatrix1, maxcorMatrix2])
maxcorMatrix = maxcorMatrix1

print maxcorMatrix

RDCs = [0, 0.2, 0.4, 0.6, 0.8, 1.]
elecRadii = [50, 100, 150, 200, 250, 300] #[400, 500, 700, 1000] #

RDCsMat, elecRadiiMat = np.meshgrid(RDCs, elecRadii)

CS = plt.contour(elecRadii, RDCs, maxcorMatrix)
plt.clabel(CS, inline=1, fontsize=10)
plt.xlabel('electrode radius [$\mu$m]')
plt.ylabel('random direction component $\\alpha$')
plt.show()

