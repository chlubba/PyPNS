# from scipy.stats import fisk
# import numpy as np
# import matplotlib.pyplot as plt
#
#
# c = 1.9
# xLogLog = np.linspace(fisk.ppf(0.01, c),fisk.ppf(0.99, c), 100)
# densityLogLog = fisk.pdf(xLogLog, c)
#
# plt.plot(xLogLog, densityLogLog)
# plt.grid()
# plt.show()

import matplotlib.pyplot as plt
import scipy
import scipy.stats
import numpy as np

# # here the distributions of myelinated and unmyelinated axon diameters are defined
# myelinatedDistribution = {
#     'densities':[100.,300.,1150.,2750.,3650.,2850.,1750.,900.,500.,250.,200.,150.,110.,100.,110.,100.,105.], #fibers densities can be given either in No/mm2 or percentage
#     'diameters': [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7.,7.5, 8., 8.5, 9.],  # corresponding diameters for each densities
# }
# unmyelinatedDistribution = {
#     'densities':[250.,1250.,5000.,8000.,9800.,10200.,8900.,7600.,5700.,4000.,3900.,2300.,2000.,1300.,900.,750.,600.,600.,500.,250.], #fibers densities can be given either in No/mm2 or percentage
#     'diameters': [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.], # corresponding diameters for each densities
# }
#
# # size = 30000
# # x = scipy.arange(size)
# # y = scipy.int(scipy.round_(scipy.stats.vonmises.rvs(5,size=size)*47))
#
# x = myelinatedDistribution['diameters']
# y = myelinatedDistribution['densities']
#
# wantedSum = 10000
# ySum = sum(y)
# scalingFactorDensity = wantedSum / ySum
# if scalingFactorDensity > 1:
#     yScaled = np.multiply(y, scalingFactorDensity)
# else:
#     yScaled = y
#
# data = np.array([], dtype=np.int64).reshape(0)
# for i in range(len(x)):
#     data = np.append(data, np.multiply(np.ones(int(yScaled[i])), x[i]))
#
# size = np.shape(data)[0]
#
# plt.plot(x,np.divide(yScaled, np.sum(yScaled)) , label='original')
# bins = np.linspace(0, 10, 100)
# # plt.hist(data, bins)
# print np.trapz(np.divide(yScaled, np.sum(yScaled)), x=x)
#
# dist_names = ['fisk', 'gamma', 'beta', 'rayleigh', 'norm', 'pareto'] # ['norm'] #
#
# for dist_name in dist_names:
#     dist = getattr(scipy.stats, dist_name)
#     param = dist.fit(data)
#
#     xFitted = np.linspace(dist.ppf(0.01, *param[:-2], loc=param[-2], scale=param[-1]), dist.ppf(0.99, *param[:-2], loc=param[-2], scale=param[-1]), 100)
#
#     pdf_fitted = dist.pdf(xFitted, *param[:-2], loc=param[-2], scale=param[-1])# * size
#
#     print np.trapz(pdf_fitted, x=xFitted)
#     plt.plot(xFitted, pdf_fitted, label=dist_name)
#     # plt.xlim(0,47)
# plt.legend(loc='upper right')
# plt.show()

# hist = np.histogram(values, nbins)
# dHistOccurence = np.double(hist[0])
# histNorm = np.divide(dHistOccurence, sum(dHistOccurence))
#
# # find the percentiles 0.01 and 0.99
# cumHistNorm = np.cumsum(histNorm)
# ip001 = np.argmax(cumHistNorm > 0.01)
# ip099 = np.argmax(cumHistNorm > 0.99)
#
#
# lengthXArray = ip099 - ip001
# p001 = hist[1][ip001]
# p099 = hist[1][ip099]
#
# xArray = np.linspace(p001, p099, lengthXArray)
#
# plt.plot(xArray, histNorm[ip001:ip099])

# here the distributions of myelinated and unmyelinated axon diameters are defined
myelinatedParams = {
    'densities':[100,300,1150,2750,3650,2850,1750,900,500,250,200,150,110,100,110,100,105], #fibers densities can be given either in No/mm2 or percentage
    'diameters': [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7.,7.5, 8., 8.5, 9.],  # corresponding diameters for each densities
}
# unkoscher correction
myelinatedParams['diameters'] = np.add(myelinatedParams['diameters'], 2.8)

myelinatedDist = {'distName' : 'manual',
                  'params' : myelinatedParams}

# myelinatedDist = {'distName' : 'normal',
#                   'params' : (5,5)}

unmyelinatedDistribution = {
    'densities':[250,1250,5000,8000,9800,10200,8900,7600,5700,4000,3900,2300,2000,1300,900,750,600,600,500,250], #fibers densities can be given either in No/mm2 or percentage
    'diameters': [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.], # corresponding diameters for each densities
}

def drawSample(distName, params, size=1):

    if distName=='constant':
        # take fixed diameter

        diam = np.multiply(np.ones(size),params)

    elif distName=='manual':
        # draw from distribution provided as diameter and density arrays

        # get diameter distribution
        densities = params['densities']
        diameters = params['diameters']

        # normalize it
        sum_d = float(sum(densities))
        normalize_densities = [x / sum_d for x in densities]

        # draw one diameter value from the distribution
        diamIndex = np.random.choice(len(normalize_densities), size, p = normalize_densities)
        diam = diameters[diamIndex]

    else:
        # draw from theoretical distribution

        dist = getattr(np.random, distName)
        diam = dist(*params, size=size)

        diam = max(0, diam)

    return diam

def getDiam(axonType):

    if axonType == 'm':

        distName = myelinatedDist['distName']
        params = myelinatedDist['params']

        drawnDiam = drawSample(distName, params, size=1)

        # choose the closest from existing axonD
        fiberD_choices = [5.7, 7.3, 8.7, 10.0, 11.5, 12.8, 14.0, 15.0, 16.0]
        diff_axonD = [abs(x - drawnDiam) for x in fiberD_choices]

        diam = fiberD_choices[np.argmin(diff_axonD)]

    elif axonType == 'u':
        distName = unmyelinatedDist['distName']
        params = unmyelinatedDist['params']

        diam = drawSample(distName, params, size=1)
    else:
        print 'Invalid axon type given to getDiam function.'
        diam = []
        quit()

    return diam

def plotDist(distName, params):

    values = drawSample(distName, params, size=1000)

    # number of bins to plot
    nbins = int(max(np.ceil(np.shape(values)[0]/10), 10))

    bins = np.linspace(min(values), max(values), nbins)
    plt.hist(values, bins)

    plt.show()


# print drawSample('lognormal', (1,1))

# plotDist('normal', (1,1))

print getDiam('m')