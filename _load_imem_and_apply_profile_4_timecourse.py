import numpy as np
import matplotlib.pylab as plt
import cPickle as pickle
import os
import matplotlib.cm as cm
import matplotlib.colors as colors


saveDict = pickle.load(open(os.path.join('/media/carl/4ECC-1C44/PyPN/imem', 'imem.dict'), "rb" ))

diameters = saveDict['diameters']
imems = saveDict['imems']
t = saveDict['t']

for i in [2]: # range(len(diameters)):

    diameter = diameters[i]
    imem = imems[i]

    # # find first and last non-zero entry
    # overX0 = [n for n, j in enumerate(imem) if j > 0.1*np.max(imem)][0]
    # overX1 = [n for n, j in enumerate(imem) if j < 0.1 * np.min(imem)][-1]
    # print t[overX1] - t[overX0]
    # # (startInd, endInd) = (a[0], a[-1])

    # alternatively, find min and max index
    startInd = np.argmax(imem)
    endInd = np.argmin(imem)

    membraneSpikeDuration = t[endInd] - t[startInd] #
    print str(membraneSpikeDuration) + ' ms'

    CV = np.sqrt(diameter) * 2

    print 'length on fibre ' + str(CV * membraneSpikeDuration) + ' mm'

    plt.plot(t, imem/np.max(imem), label='diameter' + str(diameter) + '$\mu$m')
    plt.show()

    dist_max = 100000.
    step = 2.
    positions = np.arange(-dist_max/2, dist_max/2, step)
    dts = positions * 10**-6 / CV
    dts_index = np.floor(dts/(t[1] - t[0])*10**3) # factor 1000 because ms -> s

    # sigmas = [5., 10., 50., 100., 200., 500., 1000.]
    # sigmas = [1., 2.5, 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 200., 500., 1000.]
    sigmas = np.array([0.25, 0.5, 1, 2, 3, 4])

    # descendDists = np.logspace(-4, -2, num=10)*dist_max/2
    # descendDists = np.array((0.01, 0.1, 1))*dist_max/2
    descendDists = [[200], [300], [500]]

    amp_pp = []
    for descendDist in descendDists[i]:

        # mu = 0.
        # profile = 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (positions - mu)**2 / (2 * sigma**2) )

        # descendDist = dist_max/2
        profile = np.maximum(0, np.subtract(1,np.abs(positions/descendDist)))

        # plt.plot(profile)
        # plt.show()

        # profile = 1. / (np.abs(positions)**sigma)
        profile[profile == np.inf] = np.nanmax(profile[np.logical_not(profile == np.inf)])


        profile = profile / np.max(profile)
        # plt.plot(profile, label='$\sigma$ = ' + str(sigma) + '$\mu$m')

        # plt.plot(profile)
        # plt.show()

        added_signals = np.zeros(np.shape(imem))
        colorCounter = 0
        maxAmp = 0
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

            if np.abs(position) < descendDist and np.mod(posInd,int(descendDist/50)) == 0:

                # coloring
                jet = plt.get_cmap('jet')
                cNorm = colors.Normalize(vmin=0, vmax=30)
                scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)
                colorVal = scalarMap.to_rgba(colorCounter)
                colorCounter += 1
                # print colorCounter, colorVal

                plt.plot(t, paddedSignal*profile[posInd], color = colorVal)

                maxAmp = max(maxAmp, np.max(paddedSignal*profile[posInd]))

        plt.plot(t, added_signals/np.max(added_signals)*maxAmp, color='red', linewidth=2)
        plt.title('diam ' + str(diameter) + ', cuff length ' + str(descendDist*2) + '$\mu$m')
        plt.xlim([4,7])
        # plt.show()

        # print sigma
        # print np.max(added_signals)
        # print np.min(added_signals)

        amp_pp.append(np.max(added_signals) - np.min(added_signals))

        # plt.plot(added_signals, label='$\sigma$ = ' + str(sigma) + '$\mu$m')

#     plt.semilogx(descendDists, amp_pp, label=str(diameter)+' $\mu$m') # 'diameter = ' + # /np.max(amp_pp)
# plt.xlabel('length of one triangle side ($\mu$m)')
# plt.ylabel('amplitude (unit nonsense but proportional to voltage)')
# plt.title('unmyelinated fibers SFAP amplitude vs catchment distance with linear decay (triangle)')
# plt.legend()
plt.show()