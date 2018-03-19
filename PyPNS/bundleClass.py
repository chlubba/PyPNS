from axonClass import *
import createGeometry

import numpy as np # for arrays managing
import time
import shutil
import copy

import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

from takeTime import *
import constants
import silencer

# from recordingMechanismFEMClass import RecordingMechanismFEM
# from recordingMechanismClass import RecordingMechanism
from extracellularMechanismClass import precomputedFEM

from scipy.signal import argrelextrema
from scipy.interpolate import interp1d

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

    def __init__(self, radius, numberOfAxons, pMyel, pUnmyel, length=0, paramsMyel=None, paramsUnmyel=None, axonCoords=None, bundleGuide=None,
                 segmentLengthAxon = 10, randomDirectionComponent = 0., tStop=30, timeRes=0.0025,
                 numberOfSavedSegments=300, saveV=True, saveI=False, saveLocation='Results'):

        """Constructor of the Bundle class, the main object in `PyPNS`.

        :param radius: radius of the nerve
        :param length: length of the nerve (not equal to length of the axons as they can take a curvy trajectory)
        :param bundleGuide: 3D-trajectory of the bundle

        :param numberOfAxons: number of axons exist within the nerve
        :param segmentLengthAxon: length of a straight axon segment (used by :meth:`generateGeometry.create_random_axon()`)
        :param randomDirectionComponent: (0-1) how curvy are the axons (used by :meth:`generateGeometry.create_random_axon()`)

        :param pMyel: fraction of myelinated axons
        :param pUnmyel: fraction of unmyelinated axons
        :param paramsMyel: parameters for myelinated axons
        :param paramsUnmyel: parameters for unmyelinated axons

        :param tStop: length of the simulation
        :param timeRes: either 'variable' or numerical in ms

        :param numberOfSavedSegments: specifies for how many segments the membrane voltage is saved to disk
        :param saveV: decides whether membrane voltage is saved to disk at all
        :param saveI: decides whether membrane current is saved to disk
        :param saveLocation: where to save the bundle and every other output of the calculations

        """

        self.saveI = saveI
        self.saveV = saveV

        self.paramsMyel = paramsMyel
        self.paramsUnmyel = paramsUnmyel

        self.length = length
        self.randomDirectionComponent = randomDirectionComponent
        self.segmentLengthAxon = segmentLengthAxon

        if bundleGuide is None:
            if length == 0:
                raise ValueError('Length has to be given, alternatively bundle guide.')
            self.bundleCoords = createGeometry.get_bundle_guide_straight_radius(length, segmentLengthAxon, radius=radius)
        elif bundleGuide.shape[1] == 3:
            numCoords = bundleGuide.shape[0]
            self.bundleCoords = np.hstack([bundleGuide, np.ones([numCoords,1])*radius])
        elif bundleGuide.shape[1] == 4:
            self.bundleCoords = bundleGuide
        else:
            raise ValueError('Bundle guide not valid.')

        self.excitationMechanisms = []
        self.recordingMechanisms = []

        self.pMyel = float(pMyel)/(pMyel+pUnmyel)
        self.pUnmyel = float(pUnmyel)/(pMyel+pUnmyel)

        self.numberOfAxons = numberOfAxons
        self.numMyel = int(numberOfAxons * pMyel)
        self.numUnmyel = int(numberOfAxons - self.numMyel)

        if not self.numMyel == 0 and self.paramsMyel == None:
            raise ValueError('Properties of myelinated axons not specified but at least one myelinated axon in bundle.')
        if not self.numUnmyel == 0 and self.paramsUnmyel == None:
            raise ValueError('Properties of unmyelinated axons not specified but at least one myelinated axon in bundle.')

        self.axons = []
        self.axonColors = np.zeros([self.numberOfAxons,4])
        self.radiusBundle = radius # um

        self.axonCoords = axonCoords

        # self.voltages = []
        self.trec = None  # vector of time steps

        # params for NEURON simulation
        self.tStop = tStop # set simulation duration (ms)
        self.timeRes = timeRes # set time step (ms)

        self.generate_axon_trajectories()

        self.saveParams={'timeRes': timeRes, 'tStop': tStop, 'pMyel': pMyel,
                    'paramsMyel': paramsMyel, 'paramsUnmyel': paramsUnmyel,
                    'length': length, 'numberOfAxons' : numberOfAxons, 'saveLocation': saveLocation}
        self.numberOfSavedSegments = numberOfSavedSegments

        self.basePath = get_bundle_directory(self.saveParams, new = True)

        # create axon-specific color
        jet = plt.get_cmap('Paired')
        cNorm = colors.Normalize(vmin=0, vmax=self.numberOfAxons)#len(diameters_m)-1)#
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

        # create axons
        for i in range(self.numberOfAxons):
            print "Creating axon " + str(i)
            if i < self.numUnmyel:
                self.create_axon('u', self.axonCoords[i])
            else:
                self.create_axon('m', self.axonCoords[i])
            self.axonColors[i,:] = np.array(scalarMap.to_rgba(i))


    def generate_axon_trajectories(self):
        """
        Branches on the value given as axonCoords to the bundle, with default value None. Value can be 1D-array, then
        it is considered as starting position. Or 2D, then considered as starting positions for all axons. Oor, 3D, then
        taken as the exact trajectory the axons should take.

        Returns:

        """
        axonCoordNum = len(np.shape(self.axonCoords))

        # if the user did not specify any coordinates distribute axon start positions on disk and create random trajectories
        if axonCoordNum == 0:
            axonStartPositions = self._build_disk()

            axonCoordinates = []
            for axonInd in range(self.numberOfAxons):
                c = createGeometry.create_random_axon(self.bundleCoords, axonStartPositions[axonInd,:],
                                                  self.segmentLengthAxon,
                                                  randomDirectionComponent=self.randomDirectionComponent)
                axonCoordinates.append(c)

        # if one start position is given, mainly take the bundle guide plus an offset
        elif axonCoordNum == 1:

            # only 2D input possible as it is in the y-z-plane
            assert np.shape(self.axonCoords)[0] == 2

            axonCoordinates = []
            for axonInd in range(self.numberOfAxons):
                c = createGeometry.create_random_axon(self.bundleCoords, self.axonCoords,
                                                      self.segmentLengthAxon,
                                                      randomDirectionComponent=self.randomDirectionComponent)
                axonCoordinates.append(c)

        # if two coordinates, this assumes
        elif axonCoordNum == 2:

            # one coordinate for each axon?
            assert np.shape(self.axonCoords)[0] == self.numberOfAxons
            assert np.shape(self.axonCoords)[1] == 2

            axonCoordinates = []
            for axonInd in range(self.numberOfAxons):
                c = createGeometry.create_random_axon(self.bundleCoords, self.axonCoords[axonInd],
                                                      self.segmentLengthAxon,
                                                      randomDirectionComponent=self.randomDirectionComponent)
                axonCoordinates.append(c)

        # if there is a coordinate vector for each axon
        elif axonCoordNum == 3:

            assert np.shape(self.axonCoords)[0] == self.numberOfAxons
            assert np.shape(self.axonCoords)[1] == 3

            axonCoordinates = self.axonCoords

        else:
            raise ValueError('Format of given axon coordinates not ok.')

        self.axonCoords = axonCoordinates

    def _build_disk(self):
        """Distributes the startpoints of axons uniformly over the cross section of the bundle.

        Partly from http://blog.marmakoide.org/?p=1
        """
        n = self.numberOfAxons

        radius = self.radiusBundle * np.sqrt(np.arange(n) / float(n))

        golden_angle = np.pi * (3 - np.sqrt(5))
        theta = golden_angle * np.arange(n)

        axons_pos = np.zeros((n, 2))
        axons_pos[:,0] = np.cos(theta)
        axons_pos[:,1] = np.sin(theta)
        axons_pos *= radius.reshape((n, 1))

        return axons_pos

    def create_axon(self, axonType, axonCoords):

        """
        The properties of an axon are defined. Axon type (myelinated or unmyelinated) is chosen randomly depending on
        the probabilities set in :class:Bundle. Diameters etc. are specified according to the respective dictionary
        handed to :meth:`bundleClass.Bundle.__init__`.

        :param axonType: 'u': Unmyelinated 'm': Myelinated
        """

        # # first decide by chance whether myelinated or unmyelinated
        # axonTypeIndex = np.random.choice(2,1,p = [self.pMyel, self.pUnmyel])[0]
        # axonTypes = ['m', 'u']
        # axonType = axonTypes[axonTypeIndex]

        # then get diameter. Either drawn from distribution or constant.
        axonDiameter = self._get_diam(axonType)

        # calculate the random axon coordinates
        # axonCoords = self.axonCoords[]
        # axonCoords = createGeometry.create_random_axon(self.bundleCoords, axonPosition,
        #                                                    self.segmentLengthAxon, randomDirectionComponent=self.randomDirectionComponent)

        if axonType == 'u':
            unmyel = copy.copy(self.paramsUnmyel)
            unmyel['fiberD'] = axonDiameter
            unmyel['tStop'] = self.tStop
            unmyel['timeRes'] = self.timeRes
            unmyel['numberOfSavedSegments'] = self.numberOfSavedSegments
            unmyel['rec_v'] = self.saveV
            axonParameters = dict( {'coord': axonCoords},**unmyel)
            #axonParameters = dict( {'coord': axonPosition},**unmyel)
            self.axons.append(Unmyelinated(**axonParameters))

        elif axonType == 'm':
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

        # self.axons[-1].axonPosition = axonPosition


    def _draw_sample(self, distName, params, size=1):

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
            diam = diameters[np.squeeze(diamIndex)]

        else:
            # draw from theoretical distribution

            dist = getattr(np.random, distName)
            diam = dist(*params, size=size)

            diam = max(0.3, diam)

        return diam


    def _get_diam(self, axonType):

        if axonType == 'm':

            fiberDef = self.paramsMyel['fiberD']

            if isinstance(fiberDef, dict):

                distName = fiberDef['distName']
                params = fiberDef['params']

                drawnDiam = self._draw_sample(distName, params, size=1)
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

                diam = self._draw_sample(distName, params, size=1)

            elif isinstance(fiberDef, float) or isinstance(fiberDef, int):
                diam = fiberDef
            else:
                raise Exception('Fiber diameter definition for unmyelinated axons not valid.')

            diam = max(0.1, diam)
        else:
            raise Exception('Invalid axon type given to getDiam function.')

        return diam # round(diam,1)


    def add_recording_mechanism(self, mechanism):
        """While membrane voltage and current recording can be switched off and on directly in the :class:`bundleClass.Bundle` class, extracellular recording electrodes (or other mechansisms) need to be added to the bundle before simulating with this function.

        :param mechanism: A recording mechanism instance from `recordingMechanismClass`
        """

        self.recordingMechanisms.append(mechanism)

    
    def add_excitation_mechanism(self, mechanism):
        """To incite activity on the axons containted in the nerve, excitation mechanisms need to be added to to the bundle with this function.

        :param mechanism: An excitation mechanism instance.

        """
        # mechanism.timeRes = self.timeRes
        self.excitationMechanisms.append(mechanism)


    def simulate(self):
        """
        This method starts the actual simulation of the bundle. Each axon will be created in NEURON separately, all excitation mechanisms will be connected to it and the simulation of the particular axon will be started. After having finished, the membrane currents are used to calculate the extracellular potentials for all recording mechanisms added. All NEURON objects are deleted and the next axon is processed. When all axons are done, their contributions are added in each recording mechanism to obtain the overall compound action potential.

        """

        # print '\n======================================================'
        # print '================ Starting to simulate ================'
        # print '======================================================'

        # exectue the NEURON and LFPy calculations of all axons
        self.simulate_axons()
        print '\nAll axons simulated.'

        # once simulation is over get rid of the all Neuron objects to be able to pickle the bundle-class.
        h('forall delete_section()')
        # self.trec = None
        for excitationMechanism in self.excitationMechanisms:
            excitationMechanism.delete_neuron_objects()

        # TODO: for now trec is saved as the global time vector in the bundle. However it is only valid for the CAP if a variabletime step is used.
        if isinstance(self.timeRes, numbers.Number):
            # if the time step is constant, all axons have the same time vector
            self.trec = self.axons[0].trec
        else:
            # create regularized time vector with a time step set in constants
            self.createTimeVector()

        # compute the compound action potential by summing all axon contributions up, save it, delete variables
        with takeTime('compute CAP from single axon contributions'):
            self.compute_CAPs()
        with takeTime('save CAP data'):
            self.save_CAPs_to_file()
        self.clear_CAP_vars()

        # TODO: check if size is reduces, implement more elgantly
        for recMech in self.recordingMechanisms:
            recMech.clean_up()

    def simulate_axons(self):
        """
        Routine called by :meth:`bundleClass.Bundle.simulate` to simulate all axons.
        """

        for axonIndex in range(self.numberOfAxons):

            print "\nStarting simulation of axon " + str(axonIndex)
            tStart = time.time()

            axon = self.axons[axonIndex]

            # create the neuron object specified in the axon class object
            axon.create_neuron_object()

            # connect stimulus, spontaneous spiking, etc.
            for excitationMechanism in self.excitationMechanisms:
                excitationMechanism.connect_axon(axon)

            # setup recorder for time
            # TODO: this has changed, now every axon has trec
            axon.trec = h.Vector()
            axon.trec.record(h._ref_t)

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

            # TODO: regularize time step for current and voltage here, so you don't need to care about it later
            # TODO: first check how long the interpolation takes.

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


    def save_CAPs_to_file(self):

        # t0 = time.time()

        for recMechIndex in range(len(self.recordingMechanisms)):

            recMech = self.recordingMechanisms[recMechIndex]

            # see if recording has already been performed, then do not overwirte it (with potentially empty array)
            _, CAP = self.get_CAP_from_file(recMechIndex)
            if CAP is not None:
                continue

            recMechName = recMech.__class__.__name__

            DataOut = np.array(self.trec)
            DataOut = np.column_stack( (DataOut, np.transpose(np.array(recMech.CAP))))

            filename = get_file_name('CAP_'+recMechName+'_recMech'+str(recMechIndex), self.basePath)
            # print "Save location for CAP file: " + filename

            # np.savetxt(filename, DataOut)
            np.save(filename, DataOut)

            # now save the extracellular signals of every cell
            DataOut = np.array(self.trec)
            DataOut = np.column_stack((DataOut, np.transpose(np.array(recMech.CAP_axonwise))))

            filename = get_file_name('CAP1A_'+recMechName+'_recMech'+str(recMechIndex), self.basePath)
            # print "Save location for single axon differentiated CAP file: " + filename

            # np.savetxt(filename, DataOut)
            np.save(filename, DataOut)

    def clear_all_CAP_files(self):
        """
        Get rid of all CAP recordings written to disk. This might be useful if you save the membrane current and then calculate the CAP from it several times with different recording mechanisms using :meth:`bundleClass.Bundle.compute_CAPs_from_imem_files`.
        """

        # no need to go through recording mechanisms. In fact if bundle was run but not saved, they might not even be
        # saved in bundle but still the files exist. So simply delete all directories starting with 'CAP'
        allFolderNames = [o for o in os.listdir(self.basePath) if os.path.isdir(os.path.join(self.basePath, o))]
        CAPFolderNames = [o for o in allFolderNames if o[0:3] == 'CAP']
        for CAPFolder in CAPFolderNames:
            shutil.rmtree(os.path.join(self.basePath,CAPFolder))

    def clear_all_voltage_files(self):
        """
        Get rid of all voltage recordings written to disk.
        """

        # no need to go through recordin mechanisms. In fact if bundle was run but not saved, they might not even be
        # saved in bundle but still the files exist. So simply delete all directories starting with 'CAP'
        allFolderNames = [o for o in os.listdir(self.basePath) if os.path.isdir(os.path.join(self.basePath, o))]
        CAPFolderNames = [o for o in allFolderNames if o[0:3] == 'V']
        for CAPFolder in CAPFolderNames:
            shutil.rmtree(os.path.join(self.basePath, CAPFolder))

    def clear_all_recording_mechanisms(self):
        """
        Delete all recording mechansisms from the :class:`bundleClass.Bundle` object.
        """
        self.clear_all_CAP_files()
        self.recordingMechanisms = []

    def clear_CAP_vars(self):
        for recMech in self.recordingMechanisms:
            recMech.CAP_axonwise = []
            recMech.CAP = 0

    def save_voltage_to_file_axonwise(self, axonIndex):
        """This function is called after the simulation of an individual axon and saves the voltage over time of ``numberOfSavedSegments``; ``numpy.save`` is used instead of ``numpy.savetxt`` to save space and time.

        :param axonIndex: which axon is about to be saved
        """

        # generate axon specific file name (a little clumsy, directory
        filename = get_file_name("V"+str(axonIndex), self.basePath, directoryType='V')

        # transform NEURON voltage vector to numpy array
        voltageSingleAxon = np.transpose(np.array(self.axons[axonIndex].vreclist))

        # voltageSingleAxonDownsampled = voltageSingleAxon[range(0,voltageSingleAxon.shape[0], self.downsamplingFactor), :]

        # append the sectionlength in the first column for later processing
        # in fact this is not needed now anymore because we have single files for each axon
        numberOfSegments = np.shape(voltageSingleAxon)[1]
        numberOfSegmentsArray = np.multiply(np.ones(numberOfSegments), np.array(numberOfSegments))
        voltageSingleAxonFormatted = np.row_stack((numberOfSegmentsArray, voltageSingleAxon))
        voltageSingleAxonFormatted = np.transpose(voltageSingleAxonFormatted)

        # add the time vector as the first list
        # TODO: time vector for every axon different if time step variable, save in axon and load here
        firstLine = np.transpose(np.concatenate(([0],np.array(self.axons[axonIndex].trec))))
        dataOut = np.row_stack((firstLine, voltageSingleAxonFormatted))
        # np.savetxt(filename, dataOut)
        np.save(filename, dataOut)

    def save_imem_to_file_axonwise(self, axonIndex):
        """Like :meth:`bundleClass.Bundle.save_voltage_to_file_axonwise`, the membrane currents are saved to file. All segments are saved because they might be reused for CAP calculation.

        :param axonIndex: which axon is about to be saved
        """

        # generate axon specific file name (a little clumsy, directory
        filename = get_file_name("I" + str(axonIndex), self.basePath, directoryType='I')

        # add time
        dataOut = np.row_stack((np.array(self.axons[axonIndex].trec), self.axons[axonIndex].imem))

        # np.savetxt(filename, self.axons[axonIndex].imem)
        # np.save(filename, self.axons[axonIndex].imem)
        np.save(filename, dataOut)

    def get_imem_from_file_axonwise(self, axonIndex):

        # generate axon specific file name (a little clumsy, directory
        filename = get_file_name("I" + str(axonIndex), self.basePath, directoryType='I', newFile=False)

        try:
            # imem = np.loadtxt(filename)
            imemRaw = np.load(filename)
            t = imemRaw[0, :]
            imem = imemRaw[1:,:]
            return t, imem
        except:
            raise Exception('imem-file not found for axon ' + str(axonIndex))

    def compute_CAPs_from_imem_files(self, recMecIndices=-1):
        """
        If simulation has already been run and the membrane current has been saved, this function calculates the CAP from the saved membrane current.
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
                axon.imem = self.get_imem_from_file_axonwise(axonIndex)[1]

            if len(self.recordingMechanisms) > 0:
                with takeTime("calculate extracellular potential"):
                    recMechIndex = 0
                    for recMech in self.recordingMechanisms:

                        # see if recording has already been performed
                        _, CAP = self.get_CAP_from_file(recMechIndex)
                        recMechIndex += 1
                        if CAP is not None:
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
        """Load the CAP of the recording mechansim specified from disk.

        :param recordingMechanismIndex: Index of the desired recording mechanism.

        :return: time vector, CAP
        """

        recMechName = self.recordingMechanisms[recordingMechanismIndex].__class__.__name__

        # get the whole CAP, can be single electrode or multiple
        directory = get_directory_name("CAP_"+recMechName+'_recMech'+str(recordingMechanismIndex), self.basePath)
        try:
            # newestFile = max(glob.iglob(os.path.join(directory,'')+'*.[Dd][Aa][Tt]'), key=os.path.getctime)
            newestFile = max(glob.iglob(os.path.join(directory, '') + '*.[Nn][Pp][Yy]'), key=os.path.getctime)
        except ValueError:
            # print 'No CAP calculation has been performed yet with this set of parameters.'
            return None, None

        # CAPraw = np.transpose(np.loadtxt(newestFile))
        CAPraw = np.transpose(np.load(newestFile))
        time = CAPraw[0,:]
        CAP = np.squeeze(CAPraw[1:,:])

        return time, CAP

    def get_SFAPs_from_file(self, recordingMechanismIndex=0):
        """Load the extracellular single fiber action potentials of the recording mechansim specified from disk.

        :param recordingMechanismIndex: Index of the desired recording mechanism.

        :return: time vector, matrix of SFAPs
        """

        # get the whole CAP, can be single electrode or multiple
        recMechName = self.recordingMechanisms[recordingMechanismIndex].__class__.__name__
        directory = get_directory_name('CAP1A_' + recMechName + '_recMech' + str(recordingMechanismIndex), self.basePath)
        try:
            newestFile = max(glob.iglob(os.path.join(directory, '') + '*.[Nn][Pp][Yy]'), key=os.path.getctime)
        except ValueError:
            print 'No CAP calculation has been performed yet with this set of parameters.'
            return

        # CAPraw = np.transpose(np.loadtxt(newestFile))
        SFAPsraw = np.transpose(np.load(newestFile))
        time = SFAPsraw[0, :]
        SFAPs = SFAPsraw[1:, :]

        return time, SFAPs.T

    def get_voltage_from_file(self):

        # get the voltage
        directory = get_directory_name("V", self.basePath)

        # check if calculations have been made
        try:
            # newestFile = max(glob.iglob(os.path.join(directory,'')+'*.[Dd][Aa][Tt]'), key=os.path.getctime)
            newestFile = max(glob.iglob(os.path.join(directory, '') + '*.[Nn][Pp][Yy]'), key=os.path.getctime)
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

            # Vraw = np.loadtxt(filename)
            Vraw = np.load(filename)

            timeRec = Vraw[0,1:] # extract time vector
            segmentArray = Vraw[1:,0] # extract segment numbers for each axon (varies with diameter, lambda rule)
            V = Vraw[1:,1:] # free actual voltage signals from surrounding formatting

            voltageMatrices.append(V)

        # print 'Elapsed time to load voltage file ' + str(time.time() - timeStart) + 's'
        print ("Elapsed time to load voltage file  %.2f s." % (time.time()-timeStart))

        return timeRec, voltageMatrices

    def get_voltage_from_file_one_axon(self, axonIndex):
        """Load the voltage time series of one particular axon.

        :param axonIndex: index of desired axon

        :return: time vector, voltage signal
        """

        # get axon specific file name (of existing file)
        filename = get_file_name("V" + str(axonIndex), self.basePath, newFile=False, directoryType='V')

        # Vraw = np.loadtxt(filename)
        Vraw = np.load(filename)

        timeRec = Vraw[0, 1:]  # extract time vector
        V = Vraw[1:, 1:]  # free actual voltage signals from surrounding formatting

        return timeRec, V.T

    def createTimeVector(self):
        # TODO: make this optional depending on whether a variable time step was used or not
        minTimeStep = constants.timeResResult
        minEndTime = np.Inf
        maxStartTime = -np.Inf
        for axon in self.axons:
            t = axon.trec
            # minTimeStep = np.min([minTimeStep,  np.min(np.diff(t))])
            minEndTime = np.min([minEndTime, max(t)])
            maxStartTime = np.max([maxStartTime, min(t)])

        self.trec = np.arange(maxStartTime, minEndTime, minTimeStep)

    def compute_CAPs_from_files(self):

        for recMech in self.recordingMechanisms:

            recMech.compute_CAP_from_files()

    def compute_CAPs(self):

        for recMech in self.recordingMechanisms:

            # regularize time if variable time step was used
            if not isinstance(self.timeRes, numbers.Number):
                for axonIndex, SFAP in enumerate(recMech.CAP_axonwise):
                    f = interp1d(self.axons[axonIndex].trec, SFAP)
                    recMech.CAP_axonwise[axonIndex] = f(self.trec)

            recMech.compute_overall_CAP()


