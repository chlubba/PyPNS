import cPickle as pickle
import os
import numpy as np
import matplotlib.pyplot as plt

import matplotlib.cm as cm
import matplotlib.colors as colors

import testWaveletDenoising as w
from scipy.interpolate import interp1d

saveDict = pickle.load(open(os.path.join('SFAPs', 'SFAPsPowleyMyelAsRecordings.dict'), "rb" )) # thinnerMyelDiam # SFAPsPowleyMyelAsRecordingsNewCurr

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

electrodeDistance = 70. # mm
jitterAmp = 5 #ms
jitterDist = 0.03*electrodeDistance #
numMyel = 10
numUnmyel = 10
poles = 2
poleDistance = 1 # mm
polePolarities = [1, -1]
fieldTypes = [0, 1] # 0: homo, 1: FEM

timePlotMin = 10
timePlotMax = 35

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

diametersMyel = np.array(saveDict[stringsDiam[1]])
sigma = 0.25
mu = .7
wantedNumbersOfFibers[1] = 1/(sigma * np.sqrt(2 * np.pi)) *np.exp( - (diametersMyel - mu)**2 / (2 * sigma**2) )

wantedNumbersOfFibers[0] = np.divide(wantedNumbersOfFibers[0],  np.sum(wantedNumbersOfFibers[0]))*numUnmyel
wantedNumbersOfFibers[1] = np.divide(wantedNumbersOfFibers[1],  np.sum(wantedNumbersOfFibers[0]))*numMyel

# -------------------- plot recorded data ---------------------------------
data = np.loadtxt('/Users/carl/Dropbox/_Exchange/Project/SD_1ms_AllCurrents.txt')
denoisedVoltage = w.wden(data[:,1], level=12, threshold=1.5)

tStart = 3.0245 # 3.024 #
time = data[:,0]
tCut = (time[time>tStart] - tStart)*1000
vDenCut = denoisedVoltage[time>tStart]

tCutPlot = tCut[np.logical_and(tCut>timePlotMin, tCut<timePlotMax)]
vDenCutPlot = vDenCut[np.logical_and(tCut>timePlotMin, tCut<timePlotMax)]

from scipy import signal
f_rec, Pxx_den_rec = signal.welch(vDenCutPlot - np.mean(vDenCutPlot), 1./(time[1]-time[0]), nperseg=50)


tCAP = np.arange(0,lengthOfRecording,0.0025)

fieldStrings = ['homogeneous', 'FEM']
for fieldTypeInd in [0,1]: # fieldTypes:
    CAP = np.zeros(nRecording)

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

        for fiberInd in range(numFibers):

            currentSFAP = SFAPNoArt[:, fiberInd]

            # first find maximum position in signal
            peakInd = np.argmax(currentSFAP)

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
                    wantedPeakInd = int((electrodeDistance + poleInd * poleDistance + np.random.uniform(0, 1) * jitterDist) / CV / dt)
                    # wantedPeakInd = 10/dt


                    difference = peakInd - wantedPeakInd

                    # make peak come earlier
                    if difference > 0:
                        paddedSignal = currentSFAP[difference:]
                    else:
                        paddedSignal = np.concatenate((np.zeros(np.abs(difference)), currentSFAP))

                    if len(paddedSignal) < nRecording:
                        paddedSignal = np.concatenate((paddedSignal, np.zeros(nRecording - len(paddedSignal))))
                    else:
                        paddedSignal = paddedSignal[:nRecording]

                    CAP = np.add(CAP, polePolarities[poleInd] * paddedSignal)
                    # plt.plot(paddedSignal)
                    # plt.title(str(diameters[fiberInd]))
                    # plt.show()

    # plt.plot(tCAP, CAP/np.max(CAP), label=fieldStrings[fieldTypeInd])

    f, Pxx_den = signal.welch(CAP- np.mean(CAP), 1000. / (tCAP[1] - tCAP[0]), nperseg=1024)

    # interpolate and cut to obtain same frequency range and same spacing
    f_common = f[f < np.max(f_rec)]
    interpolator_Pxx = interp1d(f_rec, Pxx_den_rec)
    Pxx_den_rec_interp = interpolator_Pxx(f_common)

    Pxx_den_cut = Pxx_den[f < np.max(f_rec)]

    # divide spectra
    Pxx_divided = Pxx_den_rec_interp / np.max(Pxx_den_rec_interp) / (Pxx_den_cut / np.max(Pxx_den_cut))

    f_finer = np.linspace(np.min(f_common), np.max(f_common), 1000)
    interpolator_Pxx_div = interp1d(f_common, Pxx_divided)
    Pxx_div_finer = interpolator_Pxx_div(f_finer)

    if fieldTypeInd == 0:
        plt.plot(f_common, 20 * np.log10(
            np.sqrt(Pxx_den_rec_interp / np.max(Pxx_den_rec_interp))), label='recorded', color='r')
    plt.plot(f_common, 20 * np.log10(
        np.sqrt(Pxx_den_cut / np.max(Pxx_den_cut))), label=fieldStrings[fieldTypeInd])#, color='b')

    # # plt.semilogy(f_common, Pxx_den_rec_interp/np.max(Pxx_den_rec_interp)/(Pxx_den_cut / np.max(Pxx_den_cut)), label='rec/sim')
    # plt.plot(f_common, 20 * np.log10(np.sqrt(Pxx_divided/ np.max(Pxx_divided))), label='rec/sim', color='g')
    # # plt.semilogy(f_finer, Pxx_div_finer, label='rec/sim')
    plt.xlabel('frequency [Hz]')
    plt.ylabel('normalized power [dB]')

    # now generate fitting filter
    from scipy import signal

    b = signal.firwin(80, 0.5, window=('kaiser', 8))
    w, h = signal.freqz(b)


plt.legend(loc='best')

plt.show()
