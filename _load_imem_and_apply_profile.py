import numpy as np
import matplotlib.pylab as plt
import cPickle as pickle
import os

saveDict = pickle.load(open(os.path.join('/media/carl/4ECC-1C44/PyPN/imem', 'imem.dict'), "rb" ))

diameters = saveDict['diameters']
imems = saveDict['imems']
t = saveDict['t']

diameter = diameters[0]
imem = imems[0]

CV = np.sqrt(diameter) * 2

dist_max = 10000
step = 2
positions = np.arange(-dist_max/2, dist_max/2, step)
dts = positions * 10**-6 * CV
dts_index = np.floor(dts/(t[1] - t[0])*10**3) # factor 1000 because ms -> s

# sigmas = [5., 10., 50., 100., 200., 500., 1000.]
# sigmas = [1., 2.5, 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 200., 500., 1000.]
sigmas = np.array([0.25, 0.5, 1, 2, 3, 4])

amp_pp = []
for sigma in sigmas:

    mu = 0.
    profile = 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (positions - mu)**2 / (2 * sigma**2) )

    # profile = 1. / (np.abs(positions)**sigma)
    profile[profile == np.inf] = np.nanmax(profile[np.logical_not(profile == np.inf)])


    profile = profile / np.max(profile)
    # plt.plot(profile, label='$\sigma$ = ' + str(sigma) + '$\mu$m')

    # plt.plot(profile)
    # plt.show()

    added_signals = np.zeros(np.shape(imem))
    for posInd, position in enumerate(positions):

        difference = dts_index[posInd]

        if difference > 0:
            paddedSignal = imem[difference:]
        else:
            paddedSignal = np.concatenate((np.zeros(np.abs(difference)), imem))

        if len(paddedSignal) < len(imem):
            paddedSignal = np.concatenate((paddedSignal, np.zeros(len(imem) - len(paddedSignal))))
        else:
            paddedSignal = paddedSignal[:len(imem)]

        added_signals += paddedSignal*profile[posInd]

    # print sigma
    # print np.max(added_signals)
    # print np.min(added_signals)

    amp_pp.append(np.max(added_signals) - np.min(added_signals))

    # plt.plot(added_signals, label='$\sigma$ = ' + str(sigma) + '$\mu$m')

plt.plot(sigmas, amp_pp)
plt.legend()
plt.show()