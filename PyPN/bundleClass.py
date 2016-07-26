import os

# load the neuron packages needed from the path of the package
PyPNDir = os.path.realpath(__file__)
PyPNDir = os.path.dirname(PyPNDir)

from neuron import h
h('load_file("noload.hoc")')
for processorSpecificFolderName in ['x86_64', 'i686', 'powerpc']:
    neuronCommand = 'nrn_load_dll("'+os.path.join(PyPNDir,processorSpecificFolderName,'.libs','libnrnmech.so')+'")'
    h(neuronCommand)

from axonClass import *
import createGeometry

import LFPy
import numpy as np # for arrays managing
import time
import shutil
import copy

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

from takeTime import *
import silencer

# from recordingMechanismFEMClass import RecordingMechanismFEM
# from recordingMechanismClass import RecordingMechanism
from LFPMechanisms import precomputedFEM

from scipy.signal import argrelextrema

from nameSetters import *

class Bundle(object):
    # radiusBundle: Radius of the bundle (typically 0.5-1.5mm)
    # numberOfAxons: Number of axons in the bundle
    # pMyel: Percentage of myelinated fiber type A
    # pUnmyel: Percentage of unmyelinated fiber type C
    #
    #
    # myelinated: parameters for fiber type A, B
    # umyelinated:  parameters for fiber type C

    def __init__(self, radius, length, numberOfAxons, pMyel, pUnmyel, paramsMyel, paramsUnmyel, bundleGuide=None,
                 segmentLengthAxon = 10, randomDirectionComponent = 0.3, tStop=30, timeRes=0.0025, numberOfSavedSegments=300, saveV=True, saveI=False, downsamplingFactor=1):

        self.saveI = saveI
        self.saveV = saveV
        self.downsamplingFactor = downsamplingFactor

        self.paramsMyel =  paramsMyel
        self.paramsUnmyel =  paramsUnmyel

        self.length = length
        self.randomDirectionComponent = randomDirectionComponent
        self.segmentLengthAxon = segmentLengthAxon

        if bundleGuide == None:
            self.bundleCoords = createGeometry.get_bundle_guide_straight_radius(length, segmentLengthAxon, radius=radius)
        elif bundleGuide.shape[1] == 3:
            numCoords = bundleGuide.shape[0]
            self.bundleCoords = np.hstack([bundleGuide, np.ones([numCoords,1])*radius])
        elif bundleGuide.shape[1] == 4:
            self.bundleCoords = bundleGuide
        else:
            raise Exception('Bundle guide not valid.')

        self.excitationMechanisms = []
        self.recordingMechanisms = []

        self.pMyel = float(pMyel)/(pMyel+pUnmyel)
        self.pUnmyel = float(pUnmyel)/(pMyel+pUnmyel)

        self.numberOfAxons = numberOfAxons
        self.axons = []
        self.axonColors = np.zeros([self.numberOfAxons,4])
        self.radiusBundle = radius # um

        # self.voltages = []

        # params for NEURON simulation
        self.tStop = tStop # set simulation duration (ms)
        self.timeRes = timeRes # set time step (ms)

        self.build_disk(self.numberOfAxons,self.bundleCoords[0,-1])

        self.saveParams={'timeRes': timeRes, 'tStop': tStop, 'pMyel': pMyel,
                    'paramsMyel': paramsMyel, 'paramsUnmyel': paramsUnmyel,
                    'length': length, 'numberOfAxons' : numberOfAxons}
        self.numberOfSavedSegments = numberOfSavedSegments

        self.basePath = get_bundle_directory(self.saveParams, new = True)

        # create axon-specific color
        jet = plt.get_cmap('Paired')
        cNorm = colors.Normalize(vmin=0, vmax=self.numberOfAxons)#len(diameters_m)-1)#
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

        # create axons
        for i in range(self.numberOfAxons):
            print "Creating axon " + str(i)

            self.create_axon(self.axons_pos[i,:])
            self.axonColors[i,:] = np.array(scalarMap.to_rgba(i))


    def build_disk(self,numberOfAxons,radiusBundle):
        """
        Partly from http://blog.marmakoide.org/?p=1
        """
        n = numberOfAxons

        radius = radiusBundle* np.sqrt(np.arange(n) / float(n))

        golden_angle = np.pi * (3 - np.sqrt(5))
        theta = golden_angle * np.arange(n)

        self.axons_pos = np.zeros((n, 2))
        self.axons_pos[:,0] = np.cos(theta)
        self.axons_pos[:,1] = np.sin(theta)
        self.axons_pos *= radius.reshape((n, 1))


    def approx_num_segs(self):

        # def lambda_f(freq, diam, Ra, cm):
        #     return 1e5 * np.sqrt(diam / (4 * np.pi *freq * Ra * cm))
        #
        # def approx_nseg_d_lambda(unmyelinatedAxon):
        #     return int((unmyelinatedAxon.L / (unmyelinatedAxon.d_lambda * lambda_f(unmyelinatedAxon.lambda_f,
        #                                                                            unmyelinatedAxon.fiberD,
        #                                                                            unmyelinatedAxon.Ra,
        #                                                                            unmyelinatedAxon.cm)) + 0.9) / 2) * 2 + 1

        nsegs_tot = 0
        nsegs_axonwise = []
        for axon in self.axons:
            if isinstance(axon, Myelinated):
                nsegTemp = axon.axontotal
            else:
                nsegTemp = axon.get_number_of_segs()
            nsegs_tot += nsegTemp

            nsegs_axonwise.append(nsegTemp)

        return nsegs_tot, nsegs_axonwise


    def create_axon(self, axonPosition):

        # first decide by chance whether myelinated or unmyelinated
        axonTypeIndex = np.random.choice(2,1,p = [self.pMyel, self.pUnmyel])
        axonTypes = ['m', 'u']
        axonType = axonTypes[axonTypeIndex]

        # then get diameter. Either drawn from distribution or constant.
        axonDiameter = self.get_diam(axonType)

        # axonCoords = np.row_stack((np.concatenate(([0], axonPosition)), np.concatenate(([self.bundleLength], axonPosition))))

        if True:
            # calculate the random axon coordinates
            axonCoords = createGeometry.create_random_axon(self.bundleCoords, axonPosition,
                                                           self.segmentLengthAxon, randomDirectionComponent=self.randomDirectionComponent)
        else:
            axonCoords = np.column_stack(axonPosition, np.concatenate(([axonPosition[0] + self.bundleLength], axonPosition[1:2])))

        if axonTypeIndex == 1:
            unmyel = copy.copy(self.paramsUnmyel)
            unmyel['fiberD'] = axonDiameter
            unmyel['tStop'] = self.tStop
            unmyel['timeRes'] = self.timeRes
            unmyel['numberOfSavedSegments'] = self.numberOfSavedSegments
            unmyel['rec_v'] = self.saveV
            axonParameters = dict( {'coord': axonCoords},**unmyel)
            #axonParameters = dict( {'coord': axonPosition},**unmyel)
            self.axons.append(Unmyelinated(**axonParameters))


        elif axonTypeIndex == 0:
            myel = copy.copy(self.paramsMyel)
            myel['fiberD'] = axonDiameter
            myel['tStop'] = self.tStop
            myel['timeRes'] = self.timeRes
            myel['numberOfSavedSegments'] = self.numberOfSavedSegments
            myel['rec_v'] = self.saveV
            axonParameters = dict( {'coord':axonCoords},**myel)
            # axonParameters = dict( {'coord':axonPosition},**myel)
            self.axons.append(Myelinated(**axonParameters))
        else:
            "Error in the draw of the axon type!"

        self.axons[-1].axonPosition = axonPosition


    def draw_sample(self, distName, params, size=1):

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


    def get_diam(self, axonType):

        if axonType == 'm':

            fiberDef = self.paramsMyel['fiberD']

            if isinstance(fiberDef, dict):

                distName = fiberDef['distName']
                params = fiberDef['params']

                drawnDiam = self.draw_sample(distName, params, size=1)
            elif isinstance(fiberDef, float) or isinstance(fiberDef, int):
                drawnDiam = fiberDef
            else:
                raise Exception('Fiber diameter definition for myelinated axons not valid.')

            # # choose the closest from existing axonD
            # fiberD_choices = [5.7, 7.3, 8.7, 10.0, 11.5, 12.8, 14.0, 15.0, 16.0]
            # diff_axonD = [abs(x - drawnDiam) for x in fiberD_choices]
            #
            # diam = fiberD_choices[np.argmin(diff_axonD)]

            diam = max(0.2, drawnDiam)

        elif axonType == 'u':

            fiberDef = self.paramsUnmyel['fiberD']

            if isinstance(fiberDef, dict):

                if 'distName' in fiberDef.keys():
                    distName = fiberDef['distName']
                    params = fiberDef['params']
                else:
                    distName = 'manual'
                    params = fiberDef

                diam = self.draw_sample(distName, params, size=1)

            elif isinstance(fiberDef, float) or isinstance(fiberDef, int):
                diam = fiberDef
            else:
                raise Exception('Fiber diameter definition for unmyelinated axons not valid.')

            diam = max(0.1, diam)
        else:
            raise Exception('Invalid axon type given to getDiam function.')

        return round(diam,1)


    def add_recording_mechanism(self, mechanism):
        self.recordingMechanisms.append(mechanism)

        # if isinstance(mechanism, RecordingMechanism):
        #     mechanism.setup_recording_elec(self.bundleCoords, self.length)
        # elif isinstance(mechanism, RecordingMechanismFEM):
        #     mechanism.bundleGuide = self.bundleCoords
        #     mechanism.setup_recording_elec(self.bundleCoords, self.length)
        # else:
        #     raise Exception('Recording mechanism type not recognized.')



    
    def add_excitation_mechanism(self, mechanism):
        mechanism.timeRes = self.timeRes
        self.excitationMechanisms.append(mechanism)


    def simulate(self):

        # exectue the NEURON and LFPy calculations of all axons
        self.simulate_axons()
        print '\nAll axons simulated.'

        # once simulation is over get rid of the all Neuron objects to be able to pickle the bundle-class.
        h('forall delete_section()')
        # self.trec = None
        for excitationMechanism in self.excitationMechanisms:
            excitationMechanism.delete_neuron_objects()

        # compute the compound action potential by summing all axon contributions up, save it, delete variables
        with takeTime('compute CAP from single axon contributions'):
            self.compute_CAPs()
        with takeTime('save CAP data'):
            self.save_CAPs_to_file()
        self.clear_CAP_vars()

    def simulate_axons(self):

        for axonIndex in range(self.numberOfAxons):

            print "\nStarting simulation of axon " + str(axonIndex)
            tStart = time.time()

            axon = self.axons[axonIndex]

            # create the neuron object specified in the axon class object
            axon.create_neuron_object()

            # connect stimulus, spontaneous spiking, etc.
            for excitationMechanism in self.excitationMechanisms:
                excitationMechanism.connect_axon(axon)

            # # todo: remove
            #
            # from scipy import signal
            # def _rectangularStimulusSignal(stimDur=0.1, amplitude=1, frequency=10, dutyCycle=0.5,
            #                                waveform='MONOPHASIC',
            #                                timeRes=.0025, delay=0, invert=False):
            #
            #     t = np.arange(0, stimDur, timeRes)
            #
            #     if waveform == 'MONOPHASIC':
            #         stimulusSignal = amplitude * 0.5 * signal.square(2 * np.pi * frequency * t,
            #                                                          duty=dutyCycle) + amplitude * 0.5
            #     elif waveform == 'BIPHASIC':
            #         stimulusSignal = amplitude * signal.square(2 * np.pi * frequency * t, duty=dutyCycle)
            #     else:
            #         print "You didn't choose the right waveform either MONOPHASIC or BIPHASIC, it has been set to default MONOPHASIC"
            #         stimulusSignal = amplitude * 0.5 * signal.square(2 * np.pi * frequency * t,
            #                                                          duty=dutyCycle) + amplitude * 0.5
            #
            #     if invert:
            #         stimulusSignal = -stimulusSignal
            #
            #     # apply delay
            #     stimulusSignal = np.concatenate((np.zeros(delay / timeRes), stimulusSignal))
            #     t = np.arange(0, stimDur + delay, timeRes)
            #
            #     return t, stimulusSignal
            #
            # t, stimulusSignal = _rectangularStimulusSignal(timeRes=self.timeRes, delay=5)
            # tVec = h.Vector(t)
            #
            # factors = np.ones(axon.totnsegs) * (-10)
            # factors = np.cumsum(factors) * axon.L / axon.totnsegs
            #
            # # global signalVectors
            # signalVectors = []
            # for sec in axon.allseclist:
            #     for segCounter, seg in enumerate(sec):
            #         stimulusSignalSection = stimulusSignal * factors[segCounter]
            #         # plt.plot(stimulusSignalSection)
            #         # plt.show()
            #
            #         signalVectors.append(h.Vector(stimulusSignalSection))
            #         signalVectors[-1].play(seg._ref_e_extracellular, tVec)
            #
            # self.trec = h.Vector()
            # self.trec.record(h._ref_t)
            #
            # axon.set_imem_recorders()
            # axon.set_voltage_recorders()
            #
            # h.celsius = axon.temperature  # set temperature in celsius
            # h.finitialize(axon.v_init)
            # h.tstop = axon.tStop
            # h.dt = axon.timeRes
            # h.run()
            #
            # axon.calc_imem()
            #
            # # todo: end remove

            # setup recorder for time
            if axonIndex == 0:
                # record time variable
                self.trec = h.Vector()
                self.trec.record(h._ref_t)

            # here we correct the conductance of the slow potassium channel from 0.08 S/cm2 to 0.12 S/cm2 to prevent
            # multiple action potentials for thin fibers
            h('forall for (x,0) if (ismembrane("axnode")) gkbar_axnode(x) = 0.12') # .16

            with takeTime("calculate voltage and membrane current"):
                axon.simulate()

            if len(self.recordingMechanisms) > 0:
                with takeTime("calculate extracellular potential"):
                    recMechIndex = 0
                    for recMech in self.recordingMechanisms:
                        recMech.compute_single_axon_CAP(axon)
                        recMechIndex += 1
            else:
                print 'No recording mechanisms added. No CAP will be recorded.'

            if self.saveV:
                with takeTime("save membrane potential to disk"):
                    self.save_voltage_to_file_axonwise(axonIndex)

            if self.saveI:
                with takeTime('save membrane current to disk'):
                    self.save_imem_to_file_axonwise(axonIndex)

            # delete the object
            axon.delete_neuron_object()

            elapsedAxon = time.time() - tStart
            # print ("Overall processing of axon %i took %.2f s. ( %.2f %% saving.)" % (axonIndex, elapsedAxon, (elapsedSaveV + elapsedSaveLFP)/elapsedAxon*100))
            print "Overall processing of axon %i took %.2f s." % (axonIndex, elapsedAxon)


    def store_geometry(self):
        self.geometry_parameters = [self.axons[0].xstart,self.axons[0].ystart,self.axons[0].zstart,self.axons[0].xend,self.axons[0].yend,self.axons[0].zend,self.axons[0].area,self.axons[0].diam,self.axons[0].length,self.axons[0].xmid,self.axons[0].ymid,self.axons[0].zmid]



    # def save_extra_rec_one_axon_one_recmech(self, axonIndex, recMechIndex):
    #
    #     directory = get_directory_name("elec" + str(recMechIndex), self.basePath)
    #
    #     # if this is the first axon to save the extracellular recording to, create the directory
    #     if axonIndex == 0:
    #         self.recordingMechanisms[recMechIndex].savePath = directory
    #         if not os.path.exists(directory):
    #             os.makedirs(directory)
    #         else:
    #             shutil.rmtree(directory)
    #             os.makedirs(directory)
    #
    #     filename = "electrode_" + str(axonIndex) + ".dat"
    #
    #     # DataOut = np.array(electrodes.LFP[0])#,1:-1
    #     # for j in range(1,len(electrodes.LFP)):
    #     #     DataOut = np.column_stack((DataOut, np.array(electrodes.LFP[j])))#,1:-1
    #
    #     # DataOut = np.transpose(np.array(electrodes.LFP))  # ,1:-1
    #     DataOut = self.recordingMechanisms[recMechIndex].CAP_axonwise[axonIndex]
    #
    #     np.savetxt(os.path.join(directory, filename), DataOut)


    # def load_one_electrode(self, elecIndex):
    #
    #     directory = get_directory_name("elec", self.basePath)
    #     filename = "electrode_"+str(elecIndex)+".dat"
    #
    #     t0 = time.time()
    #     electrodeData = np.loadtxt(os.path.join(directory, filename), unpack=True)
    #     print ("Loaded electrode %i in %.2f s." % (elecIndex, time.time()-t0))
    #
    #     return electrodeData



    def save_CAPs_to_file(self):

        # t0 = time.time()

        for recMechIndex in range(len(self.recordingMechanisms)):

            recMech = self.recordingMechanisms[recMechIndex]

            # see if recording has already been performed, then do not overwirte it (with potentially empty array)
            _, CAP = self.get_CAP_from_file(recMechIndex)
            if not CAP == None:
                continue

            recMechName = recMech.__class__.__name__

            DataOut = np.array(self.trec)
            DataOut = np.column_stack( (DataOut, np.transpose(np.array(recMech.CAP))))

            filename = get_file_name('CAP_'+recMechName+'_recMech'+str(recMechIndex), self.basePath)
            # print "Save location for CAP file: " + filename

            np.savetxt(filename, DataOut)

            # now save the extracellular signals of every cell
            DataOut = np.array(self.trec)
            DataOut = np.column_stack((DataOut, np.transpose(np.array(recMech.CAP_axonwise))))

            filename = get_file_name('CAP1A_'+recMechName+'_recMech'+str(recMechIndex), self.basePath)
            # print "Save location for single axon differentiated CAP file: " + filename

            np.savetxt(filename, DataOut)

    def clear_all_CAP_files(self):

        for recMechIndex in range(len(self.recordingMechanisms)):

            recMech = self.recordingMechanisms[recMechIndex]

            recMechName = recMech.__class__.__name__

            # filename = get_file_name("CAP_" + recMechName + '_recMech' + str(recMechIndex), self.basePath)
            directory = get_directory_name("CAP_" + recMechName + '_recMech' + str(recMechIndex), self.basePath)
            shutil.rmtree(directory)

            directory = get_directory_name('CAP1A_' + recMechName + '_recMech' + str(recMechIndex), self.basePath)
            shutil.rmtree(directory)

    def clear_all_recording_mechanisms(self):
        self.clear_all_CAP_files()
        self.recordingMechanisms = []

    def clear_CAP_vars(self):
        for recMech in self.recordingMechanisms:
            recMech.CAP_axonwise = []
            recMech.CAP = 0

    def save_voltage_to_file_axonwise(self, axonIndex):

        # generate axon specific file name (a little clumsy, directory
        filename = get_file_name("V"+str(axonIndex), self.basePath, directoryType='V')

        # transform NEURON voltage vector to numpy array
        voltageSingleAxon = np.transpose(np.array(self.axons[axonIndex].vreclist))

        voltageSingleAxonDownsampled = voltageSingleAxon[range(0,voltageSingleAxon.shape[0], self.downsamplingFactor), :]

        # append the sectionlength in the first column for later processing
        # in fact this is not needed now anymore because we have single files for each axon
        numberOfSegments = np.shape(voltageSingleAxon)[1]
        numberOfSegmentsArray = np.multiply(np.ones(numberOfSegments), np.array(numberOfSegments))
        voltageSingleAxonFormatted = np.row_stack((numberOfSegmentsArray, voltageSingleAxon))
        voltageSingleAxonFormatted = np.transpose(voltageSingleAxonFormatted)

        # add the time vector as the first list
        firstLine = np.transpose(np.concatenate(([0],np.array(self.trec))))
        dataOut = np.row_stack((firstLine, voltageSingleAxonFormatted))
        np.savetxt(filename, dataOut)

    def save_imem_to_file_axonwise(self, axonIndex):

        # generate axon specific file name (a little clumsy, directory
        filename = get_file_name("I" + str(axonIndex), self.basePath, directoryType='I')

        np.savetxt(filename, self.axons[axonIndex].imem)

    def get_imem_from_file_axonwise(self, axonIndex):

        # generate axon specific file name (a little clumsy, directory
        filename = get_file_name("I" + str(axonIndex), self.basePath, directoryType='I', newFile=False)

        try:
            imem = np.loadtxt(filename)
            return imem
        except:
            raise Exception('imem-file not found for axon ' + str(axonIndex))

    def compute_CAPs_from_imem_files(self, recMecIndices=-1):
        """
        if simulation has already been run, this function calculates the CAP from the saved membrane current
        """
        if recMecIndices == -1:
            recMecIndices = range(len(self.recordingMechanisms))
        elif not np.all([i in range(len(self.recordingMechanisms)) for i in recMecIndices]):
            print 'selected recording mechanism indices not valid. Set to all recording mechanisms.'
            recMecIndices = range(len(self.recordingMechanisms))

        for axonIndex in range(len(self.axons)):
            axon = self.axons[axonIndex]

            print ''
            with takeTime('load current of axon ' + str(axonIndex)):
                axon.imem = self.get_imem_from_file_axonwise(axonIndex)

            if len(self.recordingMechanisms) > 0:
                with takeTime("calculate extracellular potential"):
                    recMechIndex = 0
                    for recMech in self.recordingMechanisms:

                        # see if recording has already been performed
                        _, CAP = self.get_CAP_from_file(recMechIndex)
                        recMechIndex += 1
                        if not CAP == None:
                            continue

                        recMech.compute_single_axon_CAP(axon)

            else:
                print 'Nothing to do here, no recording mechanisms have been added to the bundle.'

            axon.imem = None

        # compute the compound action potential by summing all axon contributions up, save it, delete variables
        with takeTime('compute CAP from single axon contributions'):
            self.compute_CAPs()
        with takeTime('save CAP data'):
            self.save_CAPs_to_file()
        self.clear_CAP_vars()

    def get_CAP_from_file(self, recordingMechanismIndex=0):

        recMechName = self.recordingMechanisms[recordingMechanismIndex].__class__.__name__

        # get the whole CAP, can be single electrode or multiple
        directory = get_directory_name("CAP_"+recMechName+'_recMech'+str(recordingMechanismIndex), self.basePath)
        try:
            newestFile = max(glob.iglob(os.path.join(directory,'')+'*.[Dd][Aa][Tt]'), key=os.path.getctime)
        except ValueError:
            # print 'No CAP calculation has been performed yet with this set of parameters.'
            return None, None

        CAPraw = np.transpose(np.loadtxt(newestFile))
        time = CAPraw[0,:]
        CAP = np.squeeze(CAPraw[1:,:])

        return time, CAP

    def get_voltage_from_file(self):

        # get the voltage
        directory = get_directory_name("V", self.basePath)

        # check if calculations have been made
        try:
            newestFile = max(glob.iglob(os.path.join(directory,'')+'*.[Dd][Aa][Tt]'), key=os.path.getctime)
        except ValueError:
            print 'No voltage calculation has been performed yet with this set of parameters.'
            return None, None

        # load the all voltage files
        timeStart = time.time()

        # iterate over all single axon voltage files
        voltageMatrices = []
        for axonIndex in range(0,self.numberOfAxons):

            # get axon specific file name (of existing file)
            filename = get_file_name("V"+str(axonIndex), self.basePath, newFile=False, directoryType='V')

            Vraw = np.loadtxt(filename)

            timeRec = Vraw[0,1:] # extract time vector
            segmentArray = Vraw[1:,0] # extract segment numbers for each axon (varies with diameter, lambda rule)
            V = Vraw[1:,1:] # free actual voltage signals from surrounding formatting

            voltageMatrices.append(V)

        # print 'Elapsed time to load voltage file ' + str(time.time() - timeStart) + 's'
        print ("Elapsed time to load voltage file  %.2f s." % (time.time()-timeStart))

        return timeRec, voltageMatrices

    def compute_CAPs_from_files(self):

        for recMech in self.recordingMechanisms:

            recMech.compute_CAP_from_files()

    def compute_CAPs(self):

        for recMech in self.recordingMechanisms:
            recMech.compute_overall_CAP()



    def get_filename(self, recordingType):

        directory = get_directory_name(recordingType, self.basePath)
        if not os.path.exists(directory):
            os.makedirs(directory)

        # filename = 'recording.dat'
        filename = recordingType+'.dat'

        number = 0
        filenameTemp = filename
        while os.path.isfile(os.path.join(directory, filenameTemp)):
            number += 1
            # print "Be careful this file name already exist ! We concatenate a number to the name to avoid erasing your previous file."
            filenameTemp = str(number).zfill(5) + filename

        self.filename = os.path.join(directory, filenameTemp)

        return self.filename

    def conduction_velocities(self, saveToFile=False, plot=True):

        voltageMatricesTime =  self.get_voltage_from_file()

        timeArray =  voltageMatricesTime[0]
        timeStep = timeArray[1]
        voltageMatrices = voltageMatricesTime[1]


        diameterArrayMyel = []
        meanVelocityArrayMyel = []
        medianVelocityArrayMyel = []
        stdVelocityArrayMyel = []

        diameterArrayUnmyel = []
        meanVelocityArrayUnmyel = []
        medianVelocityArrayUnmyel = []
        stdVelocityArrayUnmyel = []

        # iterate over all axons
        for i in range(1,len(self.axons)):

            axon = self.axons[i]

            axonLength = axon.L
            voltageMatrix = voltageMatrices[i]

            # for myelinated axons, only the nodes will be taken into account
            if isinstance(axon, Myelinated):

                # save diameter to plot velocity against it
                diameterArrayMyel.append(axon.fiberD)

                Nnodes = self.axons[i].axonnodes
                nodePositions = range(0,(Nnodes-1)*11,11)

                # check for maximum
                maximumTimes = []
                for j in nodePositions:
                    segmentVoltage = voltageMatrix[j]

                    maximumIndex = segmentVoltage.argmax()

                    maximumTimes.append(maximumIndex)
                    plt.plot(segmentVoltage)

                # plt.title('segment voltage at nodes')
                # plt.show()

                # as only nodes get recorded, distance is node distance.
                stepSize = axon.lengthOneCycle

            else:

                # save diameter to plot velocity against it
                diameterArrayUnmyel.append(axon.fiberD)

                # number or recorded segments.
                recordedSegmentCount = np.shape(voltageMatrix)[0]

                # find the time of the signal maxima
                maximumTimes = []
                for j in range(recordedSegmentCount):#recordedSegmentCount):

                    segmentVoltage = voltageMatrix[j]

                    # # for local maxima
                    # maxima = argrelextrema(segmentVoltage, np.greater)

                    maximumIndex = segmentVoltage.argmax()

                    maximumTimes.append(maximumIndex)


                # distance between recorded segments in um
                stepSize = axonLength / recordedSegmentCount



            # time difference between action potentials of two consecutive segments in ms
            timeDifferences = np.diff(np.array(maximumTimes))*timeStep

            # conduction  velocity estimate between two segments in m/s
            velocities = stepSize/timeDifferences/1000

            # filter beginning and end
            numVel = len(velocities)
            minIndex = int(0.3*numVel)
            maxIndex = int(0.7*numVel)

            # # mask to exclude inf values (how to they evolve? stimulation artefact?)
            # maskedVelocities = np.ma.log(velocities[minIndex:maxIndex])

            # meanVelocity = maskedVelocities.mean()
            # stdVelocity = maskedVelocities.std()

            velCut = velocities[minIndex:maxIndex]

            meanVelocity = velCut.mean()
            medianVelocity = np.median(velCut)
            stdVelocity = velCut.std()

            if isinstance(axon, Myelinated):
                meanVelocityArrayMyel.append(meanVelocity)
                medianVelocityArrayMyel.append(medianVelocity)
                stdVelocityArrayMyel.append(stdVelocity)
            else:
                meanVelocityArrayUnmyel.append(meanVelocity)
                medianVelocityArrayUnmyel.append(medianVelocity)
                stdVelocityArrayUnmyel.append(stdVelocity)

            # plt.plot(velocities)
            # plt.title('axon diameter '+str(axon.fiberD)+ ' um')
            # plt.show()

        if plot:
            f, (ax1, ax2) = plt.subplots(1,2)
            ax1.scatter(diameterArrayMyel, meanVelocityArrayMyel, label='mean')
            ax1.scatter(diameterArrayMyel, medianVelocityArrayMyel, label='median', color='red')
            ax1.errorbar(diameterArrayMyel, meanVelocityArrayMyel, yerr=stdVelocityArrayMyel, linestyle='None')
            ax1.set_title('Conduction velocity of myelinated axons')
            ax1.set_xlabel('Diameter [um]')
            ax1.set_ylabel('Conduction velocity [m/s]')
            ax1.legend()

            ax2.scatter(diameterArrayUnmyel, meanVelocityArrayUnmyel, label='mean')
            ax2.scatter(diameterArrayUnmyel, medianVelocityArrayUnmyel, label='median')
            ax2.errorbar(diameterArrayUnmyel, meanVelocityArrayUnmyel, yerr=stdVelocityArrayUnmyel, linestyle='None')
            ax2.set_title('Conduction velocity of unmyelinated axons')
            ax2.set_xlabel('Diameter [um]')
            ax2.set_ylabel('Conduction velocity [m/s]')
            ax2.legend()

            # plt.show()


        # pack for return values
        myelinatedDict = {'diams': diameterArrayMyel,
                          'meanVelocity': meanVelocityArrayMyel,
                          'medianVelocity': medianVelocityArrayMyel,
                          'stdVelocity': stdVelocityArrayMyel}

        unmyelinatedDict = {'diams': diameterArrayUnmyel,
                            'meanVelocity': meanVelocityArrayUnmyel,
                          'medianVelocity': medianVelocityArrayUnmyel,
                          'stdVelocity': stdVelocityArrayUnmyel}

        returnDict = {'myel': myelinatedDict,
                      'unmyel': unmyelinatedDict}

        if saveToFile:
            pickle.dump(returnDict,open( os.path.join(self.basePath, 'conductionVelocities.dict'), "wb" ))

        return returnDict

