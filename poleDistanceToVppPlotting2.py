import numpy as np
import matplotlib.pylab as plt
import os
import matplotlib.cm as cm
import matplotlib.colors as colors

colorMapName = 'gist_stern'

radiusBundle = 150
radiusElectrode = 200
numberOfAxons = 10
unmyelinatedDiams = np.arange(0.2, 2.5, 0.3)

minDistance = 1000000
maxDistance = 0

# create recording mechanism-specific color
numDiams = len(unmyelinatedDiams)
colorMap = plt.get_cmap(colorMapName)
cNorm = colors.Normalize(vmin=0, vmax=numDiams)#len(diameters_m)-1)#
scalarMap = cm.ScalarMappable(norm=cNorm, cmap=colorMap)
diamColors = np.array(scalarMap.to_rgba(range(numDiams)))

# now load monopolar data
filename = 'vpp_against_diameter_mono'
vppsMono = np.loadtxt(os.path.join('/media/carl/4ECC-1C44/PyPN/poleDistanceAgainstVpp/Data', filename + '.txt'))

# f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, sharey=True)

# ax1.semilogy(unmyelinatedDiams, vppsMono)
# ax1.set_title('monopolar Vpp over diameter (unmyelinated)')
# ax1.set_ylabel('Vpp [mV]')
# ax1.set_xlabel('diameter [um]')

diamIndex = 0
for unmyelinatedDiam in unmyelinatedDiams:

    saveString = 'Radius Bundle %i um, radius 2D electrode %i um, %i straight unmyelinated axons of diameter %.3f um' % (
    radiusBundle, radiusElectrode, numberOfAxons, unmyelinatedDiam)

    vpps = np.loadtxt(os.path.join('/media/carl/4ECC-1C44/PyPN/poleDistanceAgainstVpp/Data', saveString+'_vpp.txt'))
    distances = np.loadtxt(os.path.join('/media/carl/4ECC-1C44/PyPN/poleDistanceAgainstVpp/Data', saveString + '_distances.txt'))

    minDistance = min(minDistance, min(distances))
    maxDistance = max(maxDistance, max(distances))

    plt.semilogy(distances, vpps, label=str(unmyelinatedDiam)+'um', color=diamColors[numDiams-1-diamIndex,:])
    plt.plot([minDistance, maxDistance], np.tile(vppsMono[diamIndex], [2,1]),  '--', color=diamColors[numDiams-1-diamIndex,:])

    # plt.show()

    diamIndex += 1

axes = plt.gca()
axes.get_xaxis().get_major_formatter().labelOnlyBase = False

plt.title('Vpp over distance of electrode poles (unmyelinated)')
plt.xlabel('electrode distance [um]')
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# #-------- myelinated axons
#
# myelinatedDiams = np.arange(0.2, 4, 0.3)
#
# minDistance = 1000000
# maxDistance = 0
#
# # now load monopolar data
# filename = 'vpp_against_diameter_mono_myelinated'
# vppsMono = np.loadtxt(os.path.join('/media/carl/4ECC-1C44/PyPN/poleDistanceAgainstVpp/Data', filename + '.txt'))
#
# ax3.semilogy(myelinatedDiams, vppsMono)
# ax3.set_title('monopolar Vpp over diameter (myelinated)')
# ax3.set_ylabel('Vpp [mV]')
# ax3.set_xlabel('diameter [um]')
#
# # create recording mechanism-specific color
# numDiams = len(myelinatedDiams)
# colorMap = plt.get_cmap(colorMapName)
# cNorm = colors.Normalize(vmin=0, vmax=numDiams)#len(diameters_m)-1)#
# scalarMap = cm.ScalarMappable(norm=cNorm, cmap=colorMap)
# diamColors = np.array(scalarMap.to_rgba(range(numDiams)))
#
# diamIndex = 0
# for myelinatedDiam in myelinatedDiams:
#     saveString = 'Radius Bundle %i um, radius 2D electrode %i um, %i straight myelinated axons of diameter %.3f um' % (
#         radiusBundle, radiusElectrode, numberOfAxons, myelinatedDiam)
#
#     vpps = np.loadtxt(os.path.join('/media/carl/4ECC-1C44/PyPN/poleDistanceAgainstVpp/Data', saveString + '_vpp.txt'))
#     distances = np.loadtxt(
#         os.path.join('/media/carl/4ECC-1C44/PyPN/poleDistanceAgainstVpp/Data', saveString + '_distances.txt'))
#
#     minDistance = min(minDistance, min(distances))
#     maxDistance = max(maxDistance, max(distances))
#
#     ax4.semilogy(distances, vpps, label=str(myelinatedDiam) + 'um', color=diamColors[numDiams-1-diamIndex, :])
#     ax4.plot([minDistance, maxDistance], np.tile(vppsMono[diamIndex], [2, 1]), '--', color=diamColors[numDiams-1-diamIndex, :])
#
#     diamIndex += 1
#
# ax4.set_title('Vpp over distance of electrode poles (myelinated)')
# ax4.set_xlabel('electrode distance [um]')
# ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.show()