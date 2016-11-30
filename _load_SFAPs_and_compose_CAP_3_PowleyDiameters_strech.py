import cPickle as pickle
import os
import numpy as np
import matplotlib.pyplot as plt

import matplotlib.cm as cm
import matplotlib.colors as colors

saveDict = pickle.load(open(os.path.join('/media/carl/4ECC-1C44/PyPN/SFAPs', 'SFAPsPowleyMyelOriginal.dict'), "rb" )) # thinnerMyelDiam


# saveDict = {'unmyelinatedDiameters' : diametersUnmyel,
#             'unmyelinatedSFAPsHomo': [],
#             'unmyelinatedSFAPsFEM': [],
#             'unmyelinatedCV' : [],
#             't': [],
#             'myelinatedDiameters': diametersMyel,
#             'myelinatedSFAPsHomo': [],
#             'myelinatedSFAPsFEM': [],
#             'myelinatedCV' : [],
#             }

# parameters
lengthOfRecording = 200 #ms
dt = 0.0025 #ms
nRecording = lengthOfRecording/dt
tArtefact = 0.1 # ms
nArtefact = tArtefact/dt

electrodeDistance = 200. # 70. # mm
jitterAmp = 5 #ms
jitterDist = 0.05*electrodeDistance # 0.03
numMyel = 1000
numUnmyel = 10000
poles = 2
poleDistance = 1 # mm
polePolarities = [1, -1]
fieldTypes = [0, 1] # 0: homo, 1: FEM

stringsDiam = ['unmyelinatedDiameters', 'myelinatedDiameters']
stringsSFAPHomo = ['unmyelinatedSFAPsHomo', 'myelinatedSFAPsHomo']
stringsSFAPFEM = ['unmyelinatedSFAPsFEM', 'myelinatedSFAPsFEM']
stringsCV = ['unmyelinatedCV', 'myelinatedCV']
tHomo = saveDict['t']
ts = [tHomo, np.arange(0,10,0.0025)]

wantedNumbersOfFibers = [(0.0691040631732923, 0.182192465406599, 0.429980837522710, 0.632957475186409, 2.05015339910575,
                          3.10696898591111,  4.54590886074274,  7.22064649366380,  7.60343269800399,  8.61543655035694,
                          8.07683524571988,  7.15617584468796,  7.04457675416097,  6.77590492234067,  5.67583310442061,
                          5.20464797635635,  3.22856301277829,  2.51011904564906,  2.06140597644239,  1.50026642131635,
                          1.32118496258518,  0.849999834520921, 0.760773515404445, 0.312027350382088, 0.200593738933586,
                          0.201222559431810, 0.15, 0.1),
                         (0.6553, 0.658, 2.245, 2.282, 4.627, 4.665, 16.734, 16.737, 19.393, 19.396, 17.776, 17.779,
                          15.503, 15.506, 9.26, 9.234, 5.993, 5.961, 2.272, 2.275, 2.138, 2.106, 1.734, 1.634, 1.151,
                          1.189, 0.948, 0.917, 2.1, 2.1)]

# diametersMyel = np.array(saveDict[stringsDiam[1]])
# sigma = 0.25
# mu = .7
# wantedNumbersOfFibers[1] =  1/(sigma * np.sqrt(2 * np.pi)) *np.exp( - (diametersMyel - mu)**2 / (2 * sigma**2) )

wantedNumbersOfFibers[0] = np.divide(wantedNumbersOfFibers[0],  np.sum(wantedNumbersOfFibers[0]))*numUnmyel
wantedNumbersOfFibers[1] = np.divide(wantedNumbersOfFibers[1],  np.sum(wantedNumbersOfFibers[0]))*numMyel

# -------------------- plot recorded data ---------------------------------
import testWaveletDenoising as w
data = np.loadtxt('/media/carl/18D40D77D40D5900/Dropbox/_Exchange/Project/SD_1ms_AllCurrents.txt')
denoisedVoltage = w.wden(data[:,1], level=12, threshold=1.5)

tStart = 3.0245 # 3.024 #
time = data[:,0]
tCut = (time[time>tStart] - tStart)*1000
vDenCut = denoisedVoltage[time>tStart]/1000

plt.plot(tCut, vDenCut, color='red', label='Experimental data')

def shift_signal(signal, difference, length):

    # make peak come earlier
    if difference > 0:
        paddedSignal = signal[difference:]
    else:
        paddedSignal = np.concatenate((np.zeros(np.abs(difference)), signal))

    if len(paddedSignal) < length:
        paddedSignal = np.concatenate((paddedSignal, np.zeros(length - len(paddedSignal))))
    else:
        paddedSignal = paddedSignal[:length]

    return paddedSignal

def shift_and_strech_signal(signal, difference, length, stretch):

    from scipy import interpolate
    f = interpolate.interp1d(range(len(signal)), signal)

    stretched_signal = f(np.arange(0,len(signal)-1,1./stretch))

    return shift_signal(stretched_signal, difference, length)

stretch = 5

tCAP = np.arange(0,lengthOfRecording,0.0025)

fieldStrings = ['Homogeneous', 'FEM']
for fieldTypeInd in [1]: # fieldTypes:

    CAP = np.zeros(nRecording)
    CAPSmoothed = np.zeros(nRecording)

    for typeInd in [1]:

        diameters = np.array(saveDict[stringsDiam[typeInd]])
        CVs = np.array(saveDict[stringsCV[typeInd]])
        t = ts[typeInd]
        numFibers = len(diameters)

        # wanted number of fibers per diameter
        # wantedNums = np.ones(numFibers)*1000
        wantedNums = wantedNumbersOfFibers[typeInd]

        if fieldTypeInd == 0:
            SFAP = np.transpose(np.array(saveDict[stringsSFAPHomo[typeInd]]))
        else:
            SFAP = np.transpose(np.array(saveDict[stringsSFAPFEM[typeInd]]))

        SFAPNoArt = SFAP [t > tArtefact, :]

        # plt.plot(SFAPNoArt)
        # plt.show()

        for fiberInd in range(numFibers):

            currentSFAP = SFAPNoArt[:, fiberInd]

            from scipy import interpolate
            f = interpolate.interp1d(range(len(currentSFAP)), currentSFAP)
            currentSFAP_streched = f(np.arange(0, len(currentSFAP) - 1, 1. / stretch))


            # caution, convolving!
            sigma = 0.6
            gx = np.arange(-3 * sigma, 3 * sigma, dt)
            gaussian = np.exp(-(gx / sigma) ** 2 / 2)
            # gaussianN = gaussian/np.sum(gaussian)
            print np.sum(gaussian)
            smoothedSFAP = np.convolve(currentSFAP, gaussian, mode="same")

            # first find maximum position in signal
            peakInd = np.argmax(currentSFAP)

            peakInd_stretched = np.argmax(currentSFAP_streched)

            # plt.plot(diameters, CVs)
            # plt.show()
            if typeInd == 1:
                CV = CVs[fiberInd]
            else:
                CV = np.sqrt(diameters[fiberInd]) * 2

            # print 'diameter : ' + str(diameters[fiberInd]) + 'CV = ' + str(CV)

            for ii in range(int(max(wantedNums[fiberInd], 1))):



                for poleInd in range(poles):

                    # wantedPeakInd = (electrodeDistance + poleInd*poleDistance) / CV / dt + jitterAmp*np.random.uniform(-1,1) / dt
                    wantedPeakInd = (electrodeDistance + poleInd * poleDistance + np.random.uniform(0, 1) * jitterDist) / CV / dt
                    # wantedPeakInd = 10/dt


                    difference = peakInd - wantedPeakInd
                    difference_stretched = peakInd_stretched - wantedPeakInd

                    paddedSignal = shift_signal(currentSFAP, difference, nRecording)
                    paddedSignalSmoothed = shift_signal(smoothedSFAP, difference, nRecording)

                    # plt.figure()
                    # plt.plot(paddedSignal)

                    paddedSignal = shift_and_strech_signal(currentSFAP, difference, nRecording, stretch=stretch)
                    paddedSignalSmoothed = shift_and_strech_signal(smoothedSFAP, difference_stretched, nRecording, stretch=stretch)

                    # plt.plot(paddedSignal)
                    # plt.show()

                    CAP = np.add(CAP, polePolarities[poleInd] * paddedSignal)
                    CAPSmoothed = np.add(CAPSmoothed, polePolarities[poleInd] * paddedSignalSmoothed)

                    # plt.plot(paddedSignal)
                    # plt.title(str(diameters[fiberInd]))
                    # plt.show()

                #     plt.plot(tCAP, CAP, label=fieldStrings[fieldTypeInd])
                # plt.show()

    # plt.plot(tCAP, CAPSmoothed, label=fieldStrings[fieldTypeInd] + ' smoothed with sigma = ' + str(sigma) + ' ms')

    scaling = np.max(vDenCut[tCut > 4])/np.max(CAP)

    plt.plot(tCAP, CAP*scaling, label=fieldStrings[fieldTypeInd] + ' scaled by ' + str(scaling))
    plt.title('Comparison between experimental data and simulation')
    plt.ylabel('$V_{ext}$ [mV]')


    # from scipy import signal
    #
    # b, a = signal.butter(4, 10 / (np.pi / dt), 'low', analog=False)
    #
    # # plt.figure()
    # # w, h = signal.freqz(b, a)
    # # plt.plot(w, 20 * np.log10(abs(h)))
    # # plt.xscale('log')
    # # plt.title('Butterworth filter frequency response')
    # # plt.xlabel('Frequency [radians / second]')
    # # plt.ylabel('Amplitude [dB]')
    # # plt.margins(0, 0.1)
    # # plt.grid(which='both', axis='both')
    # # # plt.axvline(100, color='green')  # cutoff frequency
    # # plt.show()
    #
    # y = signal.filtfilt(b, a, CAP)
    # plt.plot(tCAP, y)

    # timesTicks = np.arange(5, int(max(tCAP)), 5)
    # tickLabelStrings = []
    # for tickTime in timesTicks:
    #     tickLabelStrings.append('%2.3f' % (electrodeDistance / tickTime))
    #
    # plt.xticks(timesTicks, rotation='vertical')
    # plt.gca().set_xticklabels(tickLabelStrings)
    # plt.xlabel('conduction velocity [m/s]')
    plt.xlabel('time [ms]')

plt.grid()
plt.legend()

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