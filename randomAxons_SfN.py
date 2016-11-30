import PyPN
import matplotlib.pyplot as plt
from PyPN.takeTime import *
import numpy as np
import os

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import axes3d

# ------------------------------------------------------------------------------
# ------------------------------- SCRIPT CONTROL -------------------------------
# ------------------------------------------------------------------------------

calculationFlag = True # run simulation or load latest bundle with this parameters (not all taken into account for identification)
electricalStimulusOn = True


# ------------------------------------------------------------------------------
# --------------------------------- DEFINITION ---------------------------------
# ------------------------------------------------------------------------------

# ----------------------------- simulation params ---------------------------

tStop=20
timeRes=0.0025

# ----------------------------- bundle params -------------------------------

# set length of bundle and number of axons
lengthOfBundle = 1000 # 20000 # 400000
numberOfAxons = 30

# set the diameter distribution or fixed value
# see http://docs.scipy.org/doc/numpy/reference/routines.random.html
# 5.7, 7.3, 8.7, 10., 11.5, 12.8, 14., 15., 16.
myelinatedDiam =  0.7 # {'distName' : 'normal', 'params' : (1.0, 0.7)} # (2.0, 0.7)
unmyelinatedDiam = 5 # {'distName' : 'normal', 'params' : (0.7, 0.3)}

# axon definitions
myelinatedParameters = {'fiberD': myelinatedDiam}
unmyelinatedParameters = {'fiberD': unmyelinatedDiam}

# bundle guide
segmentLengthAxon = 30
bundleGuide = PyPN.createGeometry.get_bundle_guide_straight(lengthOfBundle, segmentLengthAxon)

# ----------------------------- stimulation params ---------------------------

# parameters of signals for stimulation
rectangularSignalParams = {'amplitude': 50.,  # Pulse amplitude (mA)
                           'frequency': 20.,  # Frequency of the pulse (kHz)
                           'dutyCycle': 0.5,  # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
                           'stimDur': 0.05,  # Stimulus duration (ms)
                           'waveform': 'MONOPHASIC',  # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
                           'delay': 2.,  # ms
                           # 'invert': True,
                           # 'timeRes': timeRes,
                           }


intraParameters = {'stimulusSignal': PyPN.signalGeneration.rectangular(**rectangularSignalParams)}

# ----------------------------- recording params -------------------------------

recordingParametersNew = {'bundleGuide': bundleGuide,
                          'radius': 200,
                          'positionAlongBundle': 7000,
                          'numberOfPoles': 2,
                          'poleDistance': 1000,
                        }

LFPMech1 = PyPN.Extracellular.precomputedFEM(bundleGuide)
LFPMech2 = PyPN.Extracellular.homogeneous(sigma=1)

electrodePos = PyPN.createGeometry.circular_electrode(**recordingParametersNew)

modularRecMech1 = PyPN.RecordingMechanism(electrodePos, LFPMech1)
modularRecMech2 = PyPN.RecordingMechanism(electrodePos, LFPMech2)

# ------------------------------------------------------------------------------
# ---------------------------------- CALCULATION -------------------------------
# ------------------------------------------------------------------------------

fig = plt.figure()

randomComponents = (0, 0.4, 0.8) # np.arange(0, 1.1, 0.2, )
numberOfRandomComponents = len(randomComponents)

for compInd, randomComponent in enumerate(randomComponents):

    # set all properties of the bundle
    bundleParameters = {'radius': 300,  # 150, #um Radius of the bundle (typically 0.5-1.5mm)
                        'length': lengthOfBundle,  # um Axon length
                        'randomDirectionComponent': randomComponent,
                        # 'bundleGuide': bundleGuide,

                        'numberOfAxons': numberOfAxons,  # Number of axons in the bundle
                        'pMyel': 0.,  # Percentage of myelinated fiber type A
                        'pUnmyel': 1.,  # Percentage of unmyelinated fiber type C
                        'paramsMyel': myelinatedParameters,  # parameters for fiber type A
                        'paramsUnmyel': unmyelinatedParameters,  # parameters for fiber type C

                        'tStop': tStop,
                        'timeRes': timeRes,

                        # 'saveI':True,
                        # 'saveV': False,

                        'numberOfSavedSegments': 50,
                        # number of segments of which the membrane potential is saved to disk
                        # 'downsamplingFactor': 100
                        }

    # create the bundle with all properties of axons and recording setup
    bundle = PyPN.Bundle(**bundleParameters)

    ax = fig.add_subplot(1, 3, compInd+1, projection='3d')

    for axonID, axon in enumerate(bundle.axons):
        ax.plot(axon.coord[:, 0], axon.coord[:, 1], axon.coord[:, 2],
                color=tuple(bundle.axonColors[axonID, :]))  # label='axon '+str(axonID),
        # ax.text(axon.coord[-1, 0], axon.coord[-1, 1], axon.coord[-1, 2], str(axonID))

    ax.plot(bundle.bundleCoords[:, 0], bundle.bundleCoords[:, 1], bundle.bundleCoords[:, 2])  # , label='bundle guide')

    plt.title(str(randomComponents[compInd]))

    # ax.xaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())
    # ax.zaxis.set_major_locator(plt.NullLocator())

    # ax.tick_params(
    #     axis='y',  # changes apply to the x-axis
    #     which='both',  # both major and minor ticks are affected
    #     bottom='off',  # ticks along the bottom edge are off
    #     top='off',  # ticks along the top edge are off
    #     labelbottom='off')  # labels along the bottom edge are off

    if False: # compInd < numberOfRandomComponents - 2:
        ax.tick_params(
            axis='x',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            bottom='off',  # ticks along the bottom edge are off
            top='off',  # ticks along the top edge are off
            labelbottom='off')  # labels along the bottom edge are off
    else:
        # ax.set_xticktabels(rotation='vertical')
        ax.set_xticks(range(0, lengthOfBundle+1001, 1000))
        ax.set_xticklabels(range(0, lengthOfBundle+1001, 1000), rotation=-30)
        # ax.set_xlabel('um')

    ax.set_zlim((-300, 300))
    ax.set_zticks((-300, 0, 300))
    ax.set_zticklabels((-300, 0, 300))#, rotation='vertical')
    # ax.set_zlabel('um')

    # ax.set_xlim(0, lengthOfBundle+1000)
    # ax.set_ylim(-(lengthOfBundle+1000), lengthOfBundle+1000)
    # ax.set_zlim(-(lengthOfBundle+1000), lengthOfBundle+1000)

    # ax.view_init(30, 269)
    ax.view_init(30, 240)

    # PyPN.plot.geometry_definition(bundle)

# plt.tight_layout()
#
# import matplotlib2tikz as mtz
# mtz.save('extrapolateMcIntyre3.tex', figureheight = '\\figureheight',
#            figurewidth = '\\figurewidth')

plt.show()

bundle = None