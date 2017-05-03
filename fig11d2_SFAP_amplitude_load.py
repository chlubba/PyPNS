import cPickle as pickle
import os
import numpy as np
import matplotlib as mpl
mpl.use('TKAgg')
import matplotlib.pyplot as plt


import matplotlib.cm as cm
import matplotlib.colors as colors

CBcolors = np.array(((0.,0.,0.), (230., 159., 0.), (86., 180., 233.), (0., 158., 115.)))/255

# saveDict = pickle.load(open(os.path.join('/media/carl/4ECC-1C44/PyPN/SFAPs', 'SFAPsOil.dict'), "rb" )) # thinnerMyelDiam
saveDict = pickle.load(open(os.path.join('SFAPs', 'SFAPsPowleyMyelAsRecordingsIdealizedCuffHigherDiams.dict'), "rb" )) # originalMyelDiam #/Volumes/SANDISK/PyPN/

amplitudes = saveDict['amplitudes']
diameters = np.arange(0.2, 4.0, 0.4)

typeStrings = ['unmyelinated', 'myelinated']
fieldStrings = ['homogeneous', 'radial inhomogeneous', 'cuff']

lineStyles = ['-', '--']

for typeInd in [0,1]:
    for fieldInd in [0,1,2]:
        plt.semilogy(diameters, amplitudes[typeInd,fieldInd,:], lineStyles[typeInd], color=CBcolors[fieldInd+1,:], label=typeStrings[typeInd]+fieldStrings[fieldInd])

plt.legend()
plt.show()