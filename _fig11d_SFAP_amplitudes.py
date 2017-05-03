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
saveDict = pickle.load(open(os.path.join('SFAPs', 'SFAPsPowleyMyelAsRecordingsIdealizedCuff2.dict'), "rb" )) # originalMyelDiam #/Volumes/SANDISK/PyPN/


from scipy.signal import butter, lfilter, freqz

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

# parameters
lengthOfRecording = 200 #ms
dt = 0.0025 #ms
nRecording = int(lengthOfRecording/dt)
tArtefact = 0.1 # ms
nArtefact = tArtefact/dt

electrodeDistance = 70.*1.2 #*1.5 # 70. # mm
jitterAmp = 5 #ms
jitterDist = 0.1*electrodeDistance # 0.03
numMyel = 100
numUnmyel = 350
poles = 2
poleDistance = 3 # 1 # mm
polePolarities = [1, -1]
fieldTypes = [0, 1] # 0: homo, 1: FEM

stringsDiam = ['unmyelinatedDiameters', 'myelinatedDiameters']
stringsSFAPHomo = ['unmyelinatedSFAPsHomo', 'myelinatedSFAPsHomo']
stringsSFAPFEM = ['unmyelinatedSFAPsFEM', 'myelinatedSFAPsFEM']
stringsSFAPIdeal = ['unmyelinatedSFAPIdeal', 'myelinatedSFAPIdeal']
stringsCV = ['unmyelinatedCV', 'myelinatedCV']

fieldStrings = ['Homogeneous', 'FEM', 'Ideal Cuff']
for fieldTypeInd in [2,1,0]: # fieldTypes:

    CAP = np.zeros(nRecording)
    CAPSmoothed = np.zeros(nRecording)

    for typeInd in [0,1]:

        diameters = np.array(saveDict[stringsDiam[typeInd]])
        CVs = np.array(saveDict[stringsCV[typeInd]])
        numFibers = len(diameters)

        if fieldTypeInd == 0:
            SFAP = np.transpose(np.array(saveDict[stringsSFAPHomo[typeInd]]))
        elif fieldTypeInd == 1:
            SFAP = np.transpose(np.array(saveDict[stringsSFAPFEM[typeInd]]))
        else:
            SFAP = np.transpose(np.array(saveDict[stringsSFAPIdeal[typeInd]]))

        plt.plot(SFAP)
        plt.show()

        SFAPNoArt = SFAP [0:, :]

        for fiberInd in range(numFibers):

            currentSFAP = SFAPNoArt[:, fiberInd]

            plt.plot(currentSFAP)
            plt.show()




# plt.grid()
plt.legend()

xlimits = ([0,120], [10,25], [50,90])
figureNames = ['CAPfull.eps', 'CAPmyel.eps', 'CAPunmyel.eps']
for limInd in range(3):
    plt.xlim(xlimits[limInd])
    if limInd == 2:
        plt.ylim((-0.01, 0.01))
    # plt.ylim(ylimits[axInd])
    # axarr[axInd].axis('equal')

    # plt.axes().set_aspect('equal', 'datalim')
    plt.savefig(os.path.join('/home/carl/Dropbox/_Exchange/Project/PyPN Paper/PythonFigureOutput', figureNames[limInd]),
            format='eps', dpi=300)

plt.show()

# jet = plt.get_cmap('jet')
# cNorm = colors.Normalize(vmin=0, vmax=numFibers - 1)
# scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)
#
# for fiberInd in range(numFibers):
#     colorVal = scalarMap.to_rgba(fiberInd)
#
#     plt.plot(t[t>tArtefact], SFAP[t>tArtefact, fiberInd], color=colorVal)
# plt.show()