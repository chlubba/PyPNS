import cPickle as pickle
import os
import numpy as np
import matplotlib as mpl
mpl.use('TKAgg')
import matplotlib.pyplot as plt


import matplotlib.cm as cm
import matplotlib.colors as colors

# saveDict = pickle.load(open(os.path.join('/media/carl/4ECC-1C44/PyPN/SFAPs', 'SFAPsOil.dict'), "rb" )) # thinnerMyelDiam
saveDict = pickle.load(open(os.path.join('SFAPs', 'SFAPsPowleyMyelAsRecordingsIdealizedCuff2.dict'), "rb" )) # originalMyelDiam #/Volumes/SANDISK/PyPN/

SFAPkeysUnmyel = ['unmyelinatedSFAPsHomo', 'unmyelinatedSFAPsFEM', 'unmyelinatedSFAPIdeal']
SFAPkeysMyel = ['myelinatedSFAPsHomo', 'myelinatedSFAPsFEM', 'myelinatedSFAPIdeal']

# unmyelSFAPs = saveDict['unmyelinatedSFAPIdeal']
myelSFAPs = saveDict['myelinatedSFAPsFEM']

f, axarr = plt.subplots(3,2, sharex=True, sharey=True)
for ind, SFAPkey in enumerate(SFAPkeysUnmyel):
    axarr[ind, 0].plot(np.transpose(saveDict[SFAPkey]))

for ind, SFAPkey in enumerate(SFAPkeysMyel):
    axarr[ind, 1].plot(np.transpose(saveDict[SFAPkey]))

plt.show()

print 'hm'