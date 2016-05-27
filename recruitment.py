import PyPN
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import cPickle as pickle
import os
import shutil

def axonsFired(bundle):

    timeRec, voltageMatrices = bundle.get_voltage_from_file()

    # now plot
    numberOfAxons = np.shape(voltageMatrices)[0]

    spiked = np.zeros(numberOfAxons)

    for i in range(numberOfAxons):
        axonIndex = i

        voltageMatrix = np.transpose(voltageMatrices[axonIndex])

        # find out whether axon is myelinated or not
        isMyelinated = not (bundle.axons[axonIndex].name == 'unmyelinated_axon')

        axonDiameter = bundle.axons[axonIndex].fiberD
        currentNumberOfSegments = np.shape(voltageMatrix)[1]
        currentAxonLength = bundle.axons[axonIndex].L

        if not isMyelinated:
            # for j in range(currentNumberOfSegments):
            #     colorVal = scalarMap.to_rgba(int(j * currentAxonLength / currentNumberOfSegments))
            #     currentAxis.plot(timeRec, voltageMatrix[:,j], color=colorVal)
            # hm = voltageMatrix[:,j]

            maxPositions = voltageMatrix.argmax(axis=0)
            maxs = voltageMatrix.max(axis=0)

            maxPositionsBig = maxPositions[maxs > -50]

            if np.size(maxPositionsBig) > np.size(maxPositions)/10:
                spiked[axonIndex] = 1

            # plt.plot(maxPositions, maxs)
            # plt.show()
            #
            # numberOfSegments = len(maxPositions)
            #
            # firstSegment = min(3000, numberOfSegments/2)
            # reached = np.add(np.where(maxPositions[(firstSegment+1):]>maxPositions[firstSegment]), firstSegment)
            # lastSegment = max(np.squeeze(reached))
            #
            # slope, _, _, _, _ = stats.linregress(range(firstSegment, lastSegment), maxPositions[firstSegment:lastSegment])
            #
            #
            # # testIndices = np.array(np.linspace(firstSegment, lastSegment, 3))
            # # testIndices = testIndices.astype(int)
            # # differences = maxPositions[testIndices[1:-1]] - maxPositions[testIndices[0:-2]]
            # #
            # # threshold = (maxPositions[lastSegment] - maxPositions[firstSegment])/10 * 0.2
            #
            # threshold = 1 # magic!
            #
            # if slope > threshold:
            #     spiked[axonIndex] = 1
            # plt.plot(maxPositions, maxs)
            # plt.plot(range(len(maxPositions)), maxPositions)
            # plt.show()

        else:
            Nnodes = bundle.axons[axonIndex].axonnodes

            # numberOfRecordingSites = np.shape(voltageMatrix)[1]

            nodePositions = range(0,(Nnodes-1)*11,11)

            nodeDistance = bundle.axons[axonIndex].lengthOneCycle

            nodeCounter = 0
            for j in nodePositions:
                # currentAxis.plot(np.array(timeRec), np.array(voltageMatrix[:,j]), color=colorVal)
                nodeCounter += 1

    return spiked

calculationFlag = True # run simulation or load latest bundle with this parameters (not all taken into account for identification)

electricalStimulusOn = True

# set simulation params
tStop=20
timeRes=0.0025#0.0025

# set length of bundle and number of axons
lengthOfBundle = 2000
numberOfAxons = 1

# loop over stuff
# analysisName = 'unmyelinated, biphasic inverted, 50um distance 2'
# unmyelinatedDiams = np.logspace(np.log(0.1)/np.log(2), np.log(3)/np.log(2), 10,  base=2)
# amplitudes = np.logspace(np.log(0.005)/np.log(2), np.log(0.2)/np.log(2), 20, base=2)
# spikingTracker = np.zeros([len(unmyelinatedDiams), len(amplitudes)])

signalShape = 'MONOPHASIC'
inverted=False
analysisName = 'unmyelinated,' + signalShape + ' ,inverted-'+str(inverted)+', 50um distance, linear'
# unmyelinatedDiams = np.logspace(np.log(0.1)/np.log(2), np.log(3)/np.log(2), 10,  base=2)
# amplitudes = np.logspace(np.log(0.005)/np.log(2), np.log(0.2)/np.log(2), 20, base=2)
unmyelinatedDiams = np.linspace(0.1, 3, 10)
amplitudes = np.linspace(0.005, 0.2, 10)*100
spikingTracker = np.zeros([len(unmyelinatedDiams), len(amplitudes)])

saveFolder = '/media/carl/4ECC-1C44/PyPN/recruitment'
fileName = analysisName + '.dat'
savePath = os.path.join(saveFolder, fileName)

diameterIndex = 0
for unmyelinatedDiam in unmyelinatedDiams:
    amplitudeIndex = 0
    for amplitude in amplitudes:

        # set the diameter distribution or fixed value
        # see http://docs.scipy.org/doc/numpy/reference/routines.random.html
        # 5.7, 7.3, 8.7, 10., 11.5, 12.8, 14., 15., 16.
        myelinatedDiam =  {'distName' : 'uniform', 'params' : (1.5, 4)} # .2 #
        # unmyelinatedDiam = 2.3 #   {'distName' : 'uniform', 'params' : (0.1, 2)} # .2 #


        # definition of the stimulation type of the axon
        cuffParameters = {      'amplitude': amplitude, # 0.005, # 0.016,#0.2,# .0001,#1.5, #0.2, # 0.004, # 10., #  # Pulse amplitude (mA)
                                'frequency': 20., # Frequency of the pulse (kHz)
                                'dutyCycle': .5, # 0.05, # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
                                'stimDur' : 0.05, # Stimulus duration (ms)
                                'waveform': signalShape, # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
                                'radius' : 50, #um
                                'timeRes' : timeRes,
                                'delay': 5, # ms
                                'invert': inverted
        }

        # MONO, inv: 0.005
        # BI, inv: 0.015
        # MONO, non-inv: 0.025
        # BI, non-inv: 0.015

        recordingParameters = { 'radius': 200,
                                # 'numberOfElectrodes': 2,
                                'positionMax': 0.5,
                                'numberOfPoles': 2
        }

        # axon parameters
        myelinatedParameters = {'fiberD': myelinatedDiam, # um Axon diameter (5.7, 7.3, 8.7, 10.0, 11.5, 12.8, 14.0, 15.0, 16.0)
        }

        # axon parameters
        unmyelinatedParameters = {'fiberD': unmyelinatedDiam, # um Axon diameter
        }

        # set all properties of the bundle
        bundleParameters = {    'radius': 150, #150, #um Radius of the bundle (typically 0.5-1.5mm)
                                'length': lengthOfBundle, # um Axon length
                                'randomDirectionComponent' : 0,

                                'numberOfAxons': numberOfAxons, # Number of axons in the bundle
                                'pMyel': 0., # Percentage of myelinated fiber type A
                                'pUnmyel': 1., #Percentage of unmyelinated fiber type C
                                'paramsMyel': myelinatedParameters, #parameters for fiber type A
                                'paramsUnmyel': unmyelinatedParameters, #parameters for fiber type C

                                'tStop' : tStop,
                                'timeRes' : timeRes,
        }

        # # combine parameters for the bundle creation
        # Parameters = dict(bundleParameters, **recordingParameters)


        # create the bundle with all properties of axons and recording setup
        bundle = PyPN.Bundle(**bundleParameters)

        # # spiking through a single electrical stimulation
        # stimulusInstance = PyPN.Stimulus(**stimulusParameters)
        # plt.plot(stimulusInstance.t, stimulusInstance.stimulusSignal)
        # plt.title('stimulus signal without delay')
        # plt.show()
        bundle.add_excitation_mechanism(PyPN.Cuff(**cuffParameters))

        # run the simulation
        bundle.simulate()

        spikingTracker[diameterIndex, amplitudeIndex] = axonsFired(bundle)[0]

        # PyPN.plot.voltage(bundle)
        # plt.show()

        shutil.rmtree(bundle.basePath)
        bundle = None

        amplitudeIndex  += 1
    diameterIndex += 1


analysisDict = {'name': analysisName,
                'diams': unmyelinatedDiams,
                'amplitudes': amplitudes,
                'spikingTracker': spikingTracker}

pickle.dump(analysisDict, open(savePath, 'wb'))