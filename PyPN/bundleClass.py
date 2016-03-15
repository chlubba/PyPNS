import os

# load the neuron packages needed from the path of the package
PyPNDir = os.path.realpath(__file__)
PyPNDir = os.path.dirname(PyPNDir)

from neuron import h
h('load_file("noload.hoc")')
for processorSpecificFolderName in ['x86_64', 'i386']:
    neuronCommand = 'nrn_load_dll("'+os.path.join(PyPNDir,processorSpecificFolderName,'.libs','libnrnmech.so')+'")'
    h(neuronCommand)

from axonClass import *
# from stimulusClass import *
# from ExcitationMechanism import *
import createGeometry

# from neuron import h
import LFPy
import numpy as np # for arrays managing
import time
import shutil
import copy

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

import silencer

from nameSetters import *

class Bundle(object):
    # radiusBundle: Radius of the bundle (typically 0.5-1.5mm)
    # numberOfAxons: Number of axons in the bundle
    # p_A: Percentage of myelinated fiber type A
    # p_C: Percentage of unmyelinated fiber type C
    #
    # numberContactPoints: Number of points on the circle constituing the cuff electrode
    # recording_elec: Position of the recording electrode along axon in um, in "BIPOLAR" case the position along axons should be given as a couple [x1,x2]
    #
    # stim_type: Stimulation type either "INTRA" or "EXTRA" for INTRA/EXTRA_cellular stimulation
    # amplitude: pulse amplitude (nA)
    # freq: frequency of the sin pulse (kHz)
    # duty_cycle: Percentage stimulus is ON for one cycl
    # stim_dur : stimulus duration (ms)
    # dur: simulation duration (ms)
    # waveform: Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
    # number_of_elecs: number of electrodes along the bundle
    # recording_type: "BIPOLAR" or "MONOPOLAR"
    #
    # myelinated_A: parameters for fiber type A
    # umyelinated:  parameters for fiber type C

    def __init__(self, radiusBundle, lengthOfBundle, bundleGuide, numberOfAxons, p_A, p_C, numberContactPoints,
                 recordingElecPos, numberElecs, myelinated_A, unmyelinated, segmentLengthAxon = 10, randomDirectionComponent = 0.3, tStop=30, timeRes=0.0025):

        self.myelinated_A =  myelinated_A
        self.unmyelinated =  unmyelinated

        self.bundleLength = lengthOfBundle
        self.randomDirectionComponent = randomDirectionComponent
        self.segmentLengthAxon = segmentLengthAxon
        self.bundleCoords = bundleGuide

        self.excitationMechanisms = []

        self.p_A = float(p_A)/(p_A+p_C)
        self.p_C = float(p_C)/(p_A+p_C)

        self.numberOfAxons = numberOfAxons
        self.axons = []
        self.axonColors = np.zeros([self.numberOfAxons,4])
        self.radiusBundle = radiusBundle # um
        self.electrodes = []
        self.voltages = []
        self.numberContactPoints = numberContactPoints
        self.numberElecs = numberElecs
        self.recordingElecPos = recordingElecPos #um

        # params for NEURON simulation
        self.tStop = tStop # set simulation duration (ms)
        self.timeRes = timeRes # set time step (ms)

        self.build_disk(self.numberOfAxons,self.radiusBundle)

        self.saveParams={'elecCount': len(self.recordingElecPos), 'dt': timeRes, 'tStop': tStop, 'p_A': self.p_A,
                    'myelinatedDiam': self.myelinated_A['fiberD'], 'unmyelinatedDiam': self.unmyelinated['fiberD'],
                    'L': self.bundleLength}

        self.basePath = get_bundle_directory(new = True, **self.saveParams)

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

    def create_axon(self, axonPosition):

        # first decide by chance whether myelinated or unmyelinated
        axonTypeIndex = np.random.choice(2,1,p = [self.p_A, self.p_C])
        axonTypes = ['m', 'u']
        axonType = axonTypes[axonTypeIndex]

        # then get diameter. Either drawn from distribution or constant.
        axonDiameter = self.get_diam(axonType)

        # axonCoords = np.row_stack((np.concatenate(([0], axonPosition)), np.concatenate(([self.bundleLength], axonPosition))))

        if True:
            # calculate the random axon coordinates
            axonCoords = createGeometry.create_random_axon(self.bundleCoords, self.radiusBundle, axonPosition,
                                                           self.segmentLengthAxon, randomDirectionComponent=self.randomDirectionComponent)
        else:
            axonCoords = np.column_stack(axonPosition, np.concatenate(([axonPosition[0] + self.bundleLength], axonPosition[1:2])))

        if axonTypeIndex == 1:
            unmyel = copy.copy(self.unmyelinated)
            unmyel['fiberD'] = axonDiameter
            unmyel['tStop'] = self.tStop
            unmyel['timeRes'] = self.timeRes
            axonParameters = dict( {'coord': axonCoords},**unmyel)
            #axonParameters = dict( {'coord': axonPosition},**unmyel)
            self.axons.append(Unmyelinated(**axonParameters))


        elif axonTypeIndex == 0:
            myel = copy.copy(self.myelinated_A)
            myel['fiberD'] = axonDiameter
            myel['tStop'] = self.tStop
            myel['timeRes'] = self.timeRes
            axonParameters = dict( {'coord':axonCoords},**myel)
            # axonParameters = dict( {'coord':axonPosition},**myel)
            self.axons.append(Myelinated(**axonParameters))
        else:
            "Error in the draw of the axon type!"

        self.axons[-1].axonPosition = axonPosition

        # self.stim = Stimulus(self.stim_type,self.axons[i], delay[i],self.stim_dur,self.amp, self.freq,self.duty_cycle, self.stim_coord, self.waveform)

    def get_diam(self, axonType):

        if axonType == 'm':
            givenDiameter = self.myelinated_A['fiberD']

            if isinstance(givenDiameter, float) or isinstance(givenDiameter, int):
                return givenDiameter
            elif isinstance(givenDiameter, dict):
                # get diameter distribution
                fiberD = self.myelinated_A['fiberD']
                densities = fiberD['densities']

                # normalize it
                sum_d = float(sum(densities))
                normalize_densities = [x / sum_d for x in densities]

                # draw one diameter value from the distribution
                draw_diam = np.random.choice(len(normalize_densities),1,p = normalize_densities)

                # why add 2.8?
                axonD = fiberD['diameters'][draw_diam]+2.8

                # choose the closest from existing axonD
                axonD_choices = [3.4,4.6,6.9,8.1,9.2,10.4,11.5,12.7] # assuming the user consider usually in the axon diameter
                diff_axonD = [abs(x-axonD) for x in axonD_choices]
                fiberD_choices = [5.7, 7.3, 8.7, 10.0, 11.5, 12.8, 14.0, 15.0, 16.0]
                fiberD = fiberD_choices[np.argmin(diff_axonD)]
                diam = fiberD
            else:
                raise('Something is wrong with your given axon diameter for myelinated axons.')

        elif  axonType == 'u':
            givenDiameter = self.unmyelinated['fiberD']

            if isinstance(givenDiameter, float) or isinstance(givenDiameter, int):
                return givenDiameter
            elif isinstance(givenDiameter, dict):
                fiberD = self.unmyelinated['fiberD']
                densities = fiberD['densities']

                sum_d = float(sum(densities))
                normalize_densities = [x / sum_d for x in densities]

                draw_diam = np.random.choice(len(normalize_densities),1,p = normalize_densities)
                D = fiberD['diameters'][draw_diam]
                diam = D
            else:
                raise('Something is wrong with your given axon diameter for unmyelinated axons.')


        else:
            raise('Wrong axon type given to function get_diam. Valid ones: u or m')

        return diam
    
    def add_excitation_mechanism(self, mechanism):
        self.excitationMechanisms.append(mechanism)
    
    def simulate(self):

        self.simulate_axons()

        # compute the compound action potential by summing all axon contributions up
        self.compute_CAP_fromfiles()
        self.save_CAP_to_file()
        self.clear_CAP_vars()

        # get rid of the all Neuron objects to be able to pickle the bundle-class.
        h('forall delete_section()')
        self.trec = None
        self.voltages = None
        self.sum_CAP = None
        for excitationMechanism in self.excitationMechanisms:
            excitationMechanism.delete_neuron_objects()

    def simulate_axons(self):

        # where are the electrodes
        # [X,Y,Z,N] = self.setup_recording_elec()
        [X,Y,Z] = self.setup_recording_elec()

        for axonIndex in range(self.numberOfAxons):

            print "\nStarting simulation of axon " + str(axonIndex)

            axon = self.axons[axonIndex]

            # where is the axon
            axonPosition = axon.axonPosition

            electrodeParameters = {         #parameters for RecExtElectrode class
                    'sigma' : 0.3,              #Extracellular potential
                    'x' : X,  #Coordinates of electrode contacts
                    'y' : Y-axonPosition[0],
                    'z' : Z-axonPosition[1],
                    # 'n' : 20,
                    # 'r' : 10,
                    # 'N' : N,
                    'method': "pointsource", #or "linesource"
                }



            # create the neuron object specified in the axon class object
            axon.create_neuron_object()

            # connect stimulus, spontaneous spiking, etc.
            for excitationMechanism in self.excitationMechanisms:
                excitationMechanism.connect_axon(axon)

            if axonIndex == 0:
            # record time variable
                self.trec = h.Vector()
                self.trec.record(h._ref_t)

            # take time of simulation
            temp = time.time()

            # actually start simulation of selected axon
            axon.simulate()

            elapsed1 = time.time()-temp
            print "Elapsed time calculate voltage and membrane current: " + str(elapsed1)

            # take time for LFPy calculation
            temp = time.time()

            # shut down the output, always errors at the end because membrane current too high
            with silencer.nostdout():
                self.electrodes.append(LFPy.recextelectrode.RecExtElectrode(axon, **electrodeParameters))

            # calculate LFP by LFPy from membrane current
            self.electrodes[axonIndex].calc_lfp()

            elapsed2 = time.time()-temp
            print "Elapsed time to calculate LFP from membrane current:" + str(elapsed2)

            self.save_electrode(axonIndex)
            self.electrodes[axonIndex]= None
            self.CAP_to_file = True

            # test if voltages can be recorded on the side
            # self.voltages.append(axon.vreclist)
            self.save_voltage_to_file_axonwise(axon.vreclist)

            # delete the object
            axon.delete_neuron_object()


    def store_geometry(self):
        self.geometry_parameters = [self.axons[0].xstart,self.axons[0].ystart,self.axons[0].zstart,self.axons[0].xend,self.axons[0].yend,self.axons[0].zend,self.axons[0].area,self.axons[0].diam,self.axons[0].length,self.axons[0].xmid,self.axons[0].ymid,self.axons[0].zmid]

    def save_electrode(self,i):
        directory = get_directory_name("elec", self.basePath)

        print "Saving extracellular potential of axon "+str(i)+" to disk."

        if i==0:
            if not os.path.exists(directory):
                os.makedirs(directory)
            else:
                shutil.rmtree(directory)
                os.makedirs(directory)
        filename = "electrode_"+str(i)+".dat"
        DataOut = np.array(self.electrodes[i].LFP[0])#,1:-1
        for j in range(1,len(self.electrodes[i].LFP)):
            DataOut = np.column_stack((DataOut, np.array(self.electrodes[i].LFP[j])))#,1:-1
        np.savetxt(directory+filename, DataOut)

    def load_one_electrode(self, elecIndex):

        directory = get_directory_name("elec", self.basePath)
        filename = "electrode_"+str(elecIndex)+".dat"

        t0 = time.time()
        electrodeData = np.loadtxt(directory + filename, unpack=True)
        print "loaded electrode "+ str(elecIndex) +  " in " + str(time.time()-t0)

        return electrodeData

    def save_CAP_to_file(self):

        DataOut = np.array(self.trec)
        DataOut = np.column_stack( (DataOut, np.transpose(np.array(self.sum_CAP))))

        # maybe add the header later. Right now we assume that the simulation is defined by the bundle object that get
        # always generated during the whole simulation. If files are opened independently from a bundle object, such a
        # header would be useful.
        # header = repr(parameters)

        filename = get_file_name("CAP", self.basePath)
        print "Save location for CAP file: " + filename

        np.savetxt(filename, DataOut)

        # now save the extracellular signals of every cell
        DataOut = np.array(self.trec)
        DataOut = np.column_stack((DataOut, np.transpose(self.AP_axonwise)))

        filename = get_file_name("CAP1A", self.basePath)
        print "Save location for single axon differentiated CAP file: " + filename

        np.savetxt(filename, DataOut)

    def clear_CAP_vars(self):
        self.AP_axonwise = None
        self.CAP = None

    def save_voltage_to_file_axonwise(self, vreclist):

        filename = get_file_name("V", self.basePath, newFile=False)

        # append voltages to file to save memory usage. Open file first with mode ab (append, binary)
        f=open(filename,'ab')

        voltageSingleAxon = np.transpose(np.array(vreclist))

        # append the sectionlength in the first column in order to differentiate different axons later
        numberOfSegments = np.shape(voltageSingleAxon)[1]
        numberOfSegmentsArray = np.multiply(np.ones(numberOfSegments), np.array(numberOfSegments))
        voltageSingleAxonFormatted = np.row_stack((numberOfSegmentsArray, voltageSingleAxon))
        voltageSingleAxonFormatted = np.transpose(voltageSingleAxonFormatted)

        if os.stat(filename).st_size == 0:
            firstLine = np.transpose(np.concatenate(([0],np.array(self.trec))))
            dataOut = np.row_stack( (firstLine, voltageSingleAxonFormatted))
            np.savetxt(f, dataOut)
        else:
            np.savetxt(f, voltageSingleAxonFormatted)

        f.close()


    def get_CAP_from_file(self):

        # get the whole CAP, can be single electrode or multiple
        directory = get_directory_name("CAP", self.basePath)
        try:
            newestFile = max(glob.iglob(directory+'*.[Dd][Aa][Tt]'), key=os.path.getctime)
        except ValueError:
            print 'No CAP calculation has been performed yet with this set of parameters.'
            quit()

        CAPraw = np.transpose(np.loadtxt(newestFile))
        time = CAPraw[0,:]
        CAP = CAPraw[1:,:]

        return time, CAP


    def get_voltage_from_file(self):
        # get the whole CAP, can be signle electrode or multiple
        directory = get_directory_name("V", self.basePath)
        try:
            newestFile = max(glob.iglob(directory+'*.[Dd][Aa][Tt]'), key=os.path.getctime)
        except ValueError:
            print 'No voltage calculation has been performed yet with this set of parameter.'
            return

        # load the raw voltage file
        timeStart = time.time()
        Vraw = np.loadtxt(newestFile)
        print 'Elapsed time to load voltage file ' + str(time.time() - timeStart) + 's'

        timeRec = Vraw[0,1:] # extract time vector
        segmentArray = Vraw[1:,0] # extract segment numbers for each axon (varies with diameter, lambda rule)
        V = Vraw[1:,1:] # free actual voltage signals from surrounding formatting

        # separate the axons, first get segment counts for each axon
        segmentNumbers = [segmentArray[0]]
        indexNextFirstEntry = segmentNumbers[0]
        while indexNextFirstEntry < len(segmentArray):
            segmentNumbers.append(segmentArray[indexNextFirstEntry])
            indexNextFirstEntry += segmentNumbers[-1]

        numberOfAxons = len(segmentNumbers)

        firstIndices = np.cumsum((np.insert(segmentNumbers,0,0)))

        voltageMatrices = []
        for i in range(numberOfAxons):
            startIndex = firstIndices[i]
            endIndex = firstIndices[i+1]-1
            voltageMatrices.append(V[startIndex:endIndex,:])

        return timeRec, voltageMatrices

    def compute_CAP_fromfiles(self):
        temp = time.time()

        print '\nCalculating CAP by summing up single axon contributions.'

        monopolar = len(self.recordingElecPos) == 1

        # variable to save the sum over all axons
        self.sum_CAP = np.zeros((self.numberElecs,len(self.trec)))

        # variable to save the extracellular signal from each cell separately, at the last electrode position.
        self.AP_axonwise = np.zeros((self.numberOfAxons, len(self.trec)))

        # load the recordings for every axon one by one and add them.
        for elecIndex in range(self.numberOfAxons):
            electrodeData = self.load_one_electrode(elecIndex)

        # The contactpoints that constitute one cuff electrode ring have to be recovered, summed up together per
        # recording location along the axon
            for i in range(self.numberElecs):
                if monopolar:
                    contactPointIndices = range(self.numberContactPoints*i, self.numberContactPoints*(1+i))
                    sumOverContactPoints = np.sum(electrodeData[contactPointIndices, :], 0)
                else:
                    contactPointIndicesPole1 = range(self.numberContactPoints*2*i, self.numberContactPoints*(1+2*i))
                    contactPointIndicesPole2 = range(self.numberContactPoints*(2*i+1), self.numberContactPoints*(2*(i+1)))
                    sumOverContactPoints = np.sum(electrodeData[contactPointIndicesPole1, :] - electrodeData[contactPointIndicesPole2, :], 0)

                self.sum_CAP[i,:] = self.sum_CAP[i,:] +  sumOverContactPoints

                if i == self.numberElecs-1:
                    self.AP_axonwise[elecIndex,:] = sumOverContactPoints

        elapsed = time.time()-temp
        print "Elapsed time to compute CAP " + str(elapsed) + " \n"


    def setup_recording_elec(self):

        if True:
            # calculte recording electrode positions for a 3D shaped bundle
            electrodePositions = createGeometry.electrode_positions_bundle_guided(self.bundleCoords, self.radiusBundle,
                                                                               self.numberElecs, self.numberContactPoints,
                                                                               self.recordingElecPos)
            X, Y, Z = electrodePositions[:,0], electrodePositions[:,1], electrodePositions[:,2]
        else:

            if (self.numberElecs == 1):
                # if one recording site
                if len(self.recordingElecPos) == 1:
                    # if monopolar
                    X = np.zeros(self.numberContactPoints)+self.recordingElecPos
                elif len(self.recordingElecPos) == 2:
                    # if bipolar
                    X = np.repeat(self.recordingElecPos,self.numberContactPoints,axis=0)
                elif len(self.recordingElecPos) > 2 or len(self.recordingElecPos) == 0:
                    # if wrong number of poles entered
                    raise TypeError("Only monopolar and bipolar recording are supported")
            else:
                if len(self.recordingElecPos) >1:
                    raise TypeError("Please use only the monopolar configuration of 'recording_elec' to record along the bundle")
                X1 = [np.linspace(0, self.recordingElecPos[0], self.numberElecs)]
                X = np.repeat(X1,self.numberContactPoints, axis=0)
                X = X.flatten()
            angles = np.linspace(0,360, self.numberContactPoints, endpoint = False)
            Y1 = np.round(self.radiusBundle*np.cos(angles*np.pi/180.0),2)
            Z1 = np.round(self.radiusBundle*np.sin(angles*np.pi/180.0),2)
            Y = np.tile(Y1,self.numberElecs*len(self.recordingElecPos))
            Z = np.tile(Z1,self.numberElecs*len(self.recordingElecPos))
            N = np.empty((self.numberContactPoints*self.numberElecs*len(self.recordingElecPos), 3))
            for i in xrange(N.shape[0]):
                N[i,] = [1, 0, 0] #normal vec. of contacts

        self.electrodeCoords = np.transpose(np.row_stack((X,Y,Z)))

        return [X,Y,Z]#,N]

    

    def get_filename(self, recordingType):

        directory = get_directory_name(recordingType, self.basePath)
        if not os.path.exists(directory):
            os.makedirs(directory)

        # filename = 'recording.dat'
        filename = recordingType+'.dat'

        number = 0
        filenameTemp = filename
        while os.path.isfile(directory+filenameTemp):
            number += 1
            # print "Be careful this file name already exist ! We concatenate a number to the name to avoid erasing your previous file."
            filenameTemp = str(number).zfill(5) + filename

        self.filename = directory+filenameTemp

        return self.filename

    