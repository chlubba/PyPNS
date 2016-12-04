import numpy as np
import matplotlib.pylab as plt
import cPickle as pickle
import os

saveDict = pickle.load(open(os.path.join('/media/carl/4ECC-1C44/PyPN/imem', 'imem.dict'), "rb" ))

diameters = saveDict['diameters']
imems = saveDict['imems']
t = saveDict['t']

for i in range(len(diameters)):

    diameter = diameters[i]
    imem = imems[i]

    # # find first and last non-zero entry
    # a = [n for n, j in enumerate(-imem) if j > 0.001] # [0,-1]
    # (startInd, endInd) = (a[0], a[-1])

    # alternatively, find min and max index
    startInd = np.argmax(imem)
    endInd = np.argmin(imem)

    membraneSpikeDuration = t[endInd] - t[startInd] #
    print str(membraneSpikeDuration) + ' ms'

    CV = np.sqrt(diameter) * 2

    print 'length on fibre ' + str(CV * membraneSpikeDuration) + ' mm'

    # plt.plot(t*CV, imem)
    # plt.show()



    dist_max = 100000.
    step = 2.
    positions = np.arange(-dist_max/2, dist_max/2, step)
    dts = positions * 10**-6 / CV
    dts_index = np.floor(dts/(t[1] - t[0])*10**3) # factor 1000 because ms -> s

    # sigmas = [5., 10., 50., 100., 200., 500., 1000.]
    # sigmas = [1., 2.5, 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 200., 500., 1000.]
    sigmas = np.array([0.25, 0.5, 1, 2, 3, 4])

    descendDists = np.logspace(-4, -0, num=20)*dist_max/2
    # descendDists = np.array((0.01, 0.1, 1))*dist_max/2

    amp_pp = []
    amp_int = []
    for descendDist in descendDists:

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
        amp_int.append(np.sum(np.abs(added_signals)))

        # plt.plot(added_signals, label='$\sigma$ = ' + str(sigma) + '$\mu$m')

    plt.loglog(descendDists*2, amp_pp, label=str(diameter)+' $\mu$m') # 'diameter = ' + # /np.max(amp_pp)
    # plt.semilogx(descendDists, amp_int, label=str(diameter) + ' $\mu$m')  # 'diameter = ' + # /np.max(amp_pp)
plt.xlabel('length of cuff/ oil ($\mu$m)')
plt.ylabel('amplitude (proportional to voltage)')
plt.title('unmyelinated fibers SFAP amplitude vs oil length')
plt.legend()
plt.show()