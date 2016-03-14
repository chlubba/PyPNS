from axonClass import *
from upstreamSpikingClass import *
from stimulusClass import *
import createGeometry

# from neuron import h
import refextelectrode
import numpy as np # for arrays managing
import math
import time
import glob
import os
import shutil
import copy

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d



from nameSetters import getDirectoryName, getFileName, getBundleDirectory

class Bundle(object):
    # radius_bundle: Radius of the bundle (typically 0.5-1.5mm)
    # number_of_axons: Number of axons in the bundle
    # p_A: Percentage of myelinated fiber type A
    # p_C: Percentage of unmyelinated fiber type C
    #
    # number_contact_points: Number of points on the circle constituing the cuff electrode
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

    def __init__(self, radius_bundle, lengthOfBundle, segmentLengthAxon, bundleGuide, number_of_axons, p_A, p_C, number_contact_points,
                 recording_elec_pos, jitter_para, stim_type, duty_cycle, freq, amplitude, stim_dur, dur, number_elecs,
                 myelinated_A, unmyelinated,rec_CAP, waveform, randomDirectionComponent = 0.3):

        self.myelinated_A =  myelinated_A
        self.unmyelinated =  unmyelinated

        self.bundleLength = lengthOfBundle
        self.randomDirectionComponent = randomDirectionComponent
        self.segmentLengthAxon = segmentLengthAxon
        self.bundleCoords = bundleGuide

        self.waveform = waveform
        self.stim_type = stim_type
        self.stim_coord = [[0, radius_bundle*math.cos(math.pi/number_contact_points),
                            radius_bundle*math.sin(math.pi/number_contact_points)]]
        self.duty_cycle = duty_cycle
        self.freq = freq
        self.amp = amplitude
        self.stim_dur = stim_dur
        self.dur = dur
        self.rec_CAP = rec_CAP

        self.p_A = float(p_A)/(p_A+p_C)
        self.p_C = float(p_C)/(p_A+p_C)

        self.number_of_axons = number_of_axons
        self.axons = []
        self.axonColors = np.zeros([self.number_of_axons,4])
        self.radius_bundle = radius_bundle # um
        self.electrodes = []
        self.voltages = []
        self.number_contact_points = number_contact_points
        self.number_elecs = number_elecs
        self.recording_elec_pos = recording_elec_pos #um

        self.build_disk(self.number_of_axons,self.radius_bundle)

        self.saveParams={'elecCount': len(self.recording_elec_pos), 'dt': h.dt, 'tStop': h.tstop, 'p_A': self.p_A,
                    'myelinatedDiam': self.myelinated_A['fiberD'], 'unmyelinatedDiam': self.unmyelinated['diam'],
                    'L': self.bundleLength, 'stimType': self.stim_type, 'stimWaveform' : self.waveform,
                    'stimDutyCycle': self.duty_cycle, 'stimAmplitude' : self.amp}

        self.basePath = getBundleDirectory(new = True, **self.saveParams)


        # maybe place this somewhere else, but for now:
        # create axon-specific color
        jet = plt.get_cmap('Paired')
        cNorm = colors.Normalize(vmin=0, vmax=self.number_of_axons)#len(diameters_m)-1)#
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)


        # create axons within one fourth slice of the whole bundle.
        self.virtual_number_axons = 0
        for i in range(self.number_of_axons):
            # if within the desired subsection of the bundle, create axon
            if True:#((self.axons_pos[i,1]>=0 and self.axons_pos[i,1]< self.axons_pos[i,0]) or (self.axons_pos[i,0]== 0 and self.axons_pos[i,1]==0)):
                print "Creating axon " + str(i)

                self.create_axon(self.axons_pos[i,:])
                self.axonColors[i,:] = np.array(scalarMap.to_rgba(i))
                self.virtual_number_axons +=1


        if not self.stim_type == 'NONE':
            # create Simulus instace used for all axons
            self.stim = Stimulus(self.stim_type, self.stim_dur,self.amp, self.freq,self.duty_cycle, self.stim_coord, self.waveform)


    def addUpstreamSpiking(self, tStart=0., tStop=h.tstop, lambd = 1000., correlation = 0.1):
        # create upstream activity
        self.upstreamSpiking = UpstreamSpiking(self.number_of_axons, tStart=tStart, tStop=tStop, lambd=lambd, correlation=correlation)

    def simulateBundle(self):

        self.simulateAxons()

        if self.rec_CAP:
            self.compute_CAP_fromfiles()
            self.save_CAP_to_file()
            self.clear_CAP_vars()

        # axonLimit = 15
        # if self.virtual_number_axons <= axonLimit:
        #     self.save_voltage_to_file()
        # else:
        #    print 'Voltage not saved, too many axons (' +str(self.virtual_number_axons) + ', but maximum ' + str(axonLimit) + ').'

        # get rid of the all Neuron objects to be able to pickle the bundle-class.
        h('forall delete_section()')
        self.trec = None
        if not self.stim_type == 'NONE':
            self.stim.svec = None
        for axon in self.axons:
            # from Stimulus
            axon.stim = None
            # from upstreamSpiking
            axon.synapse = None
            axon.netCon = None
            axon.spikeVec = None
            axon.vecStim = None

            try:
                for memirec in axon.memireclist:
                    memirec = None
                axon.memireclist = None
            except:
                pass
            try:
                for vrec in axon.vreclist:
                    vrec = None
                axon.vreclist = None
            except:
                pass
            axon.allseclist = None
            # also delete unnecessary data that will no longer be used to keep the pickled file small
            axon.imem = None
        self.voltages = None
        self.sum_CAP = None



    def save_CAP_to_file(self):

        DataOut = np.array(self.trec)
        DataOut = np.column_stack( (DataOut, np.transpose(np.array(self.sum_CAP))))

        # maybe add the header later. Right now we assume that the simulation is defined by the bundle object that get
        # always generated during the whole simulation. If files are opened independently from a bundle object, such a
        # header would be useful.
        # header = repr(parameters)

        filename = getFileName("CAP", self.basePath)
        print "Save location for CAP file: " + filename

        np.savetxt(filename, DataOut)

        # now save the extracellular signals of every cell
        DataOut = np.array(self.trec)
        DataOut = np.column_stack((DataOut, np.transpose(self.AP_axonwise)))

        filename = getFileName("CAP1A", self.basePath)
        print "Save location for single axon differentiated CAP file: " + filename

        np.savetxt(filename, DataOut)

    def clear_CAP_vars(self):
        self.AP_axonwise = None
        self.CAP = None

    def save_voltage_to_file_axonwise(self, vreclist):

        filename = getFileName("V", self.basePath, newFile=False)

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
        directory = getDirectoryName("CAP", self.basePath)
        try:
            newestFile = max(glob.iglob(directory+'*.[Dd][Aa][Tt]'), key=os.path.getctime)
        except ValueError:
            print 'No CAP calculation has been performed yet with this set of parameters.'
            quit()

        CAPraw = np.transpose(np.loadtxt(newestFile))
        time = CAPraw[0,:]
        CAP = CAPraw[1:,:]

        return time, CAP

    def plot_geometry(self):

        if len(self.axons[0].xmid) == 0:
            print 'Bundle has not been run yet. No geometry information was generated in NEURON.'
            return

        fig = plt.figure()
        ax = fig.gca(projection='3d')

        axonID = 0
        for axon in self.axons:
            if type(axon) == Myelinated:
                style = '--'
            else:
                style = '-'
            ax.plot(axon.xmid, axon.ymid, axon.zmid, style, label='axon '+str(axonID), color= tuple(self.axonColors[axonID,:]))
            ax.text(axon.xmid[-1], axon.ymid[-1], axon.zmid[-1], str(axonID))
            axonID += 1
        ax.plot(self.bundleCoords[:,0], self.bundleCoords[:,1], self.bundleCoords[:,2], label='bundle guide')
        plt.legend()

        elecCoords = self.electrodeCoords
        ax.scatter(elecCoords[:,0], elecCoords[:,1], elecCoords[:,2])

        elecPoles = len(self.recording_elec_pos)
        for i in range(self.number_elecs):
            for j in range(elecPoles):
                # selectionIndices = range(i+j*self.number_contact_points, self.number_contact_points*self.number_elecs + j*self.number_contact_points, self.number_elecs)
                selectionIndices = range((i*elecPoles+j)*self.number_contact_points, self.number_contact_points*(i*elecPoles+j+1))

                ringCoords = elecCoords[selectionIndices,:]
                ringCoords = np.row_stack((ringCoords, ringCoords[0,:]))
                ax.plot(ringCoords[:,0], ringCoords[:,1], ringCoords[:,2], color=[0.8,0.8,0.8])


        plt.savefig(self.basePath+'geometry.png')

    def plot_CAP1D_singleAxon(self, maxNumberOfAxons):

        # get the whole CAP, can be single electrode or multiple
        directory = getDirectoryName("CAP1A", self.basePath)
        try:
            newestFile = max(glob.iglob(directory+'*.[Dd][Aa][Tt]'), key=os.path.getctime)
        except ValueError:
            print 'No CAP calculation has been performed yet with this set of parameters.'
            quit()

        CAPraw = np.transpose(np.loadtxt(newestFile))
        time = CAPraw[0,:]
        CAP = CAPraw[1:,:]

        numberOfPlots = min(maxNumberOfAxons, self.virtual_number_axons)
        axonSelection = np.floor(np.linspace(0,self.virtual_number_axons-1, numberOfPlots))

        # Subplots
        f, axarr = plt.subplots(numberOfPlots, sharex=True)

        for i in range(len(axonSelection)):
            axonIndex = int(axonSelection[i])

            axon = self.axons[i]

            axonDiameter = axon.fiberD

            if type(self.axons[axonIndex]) == Myelinated:
                axonType = 'myelinated'
            else:
                axonType = 'unmyelinated'

            CAPSingleAxon = CAP[axonIndex,:]

            axarr[i].plot(time, CAPSingleAxon, color= tuple(self.axonColors[axonIndex,:]))
            axarr[i].set_title('Axon ' + str(axonIndex) + ' (' + axonType + ') with diameter ' + str(axonDiameter) + 'um')
            axarr[i].set_ylabel('CAP [mV]')

        plt.savefig(self.basePath+'CAPSingleAxons.png')





    def plot_CAP1D(self, maxNumberOfSubplots = 10):

        # first load the desired data from file
        time, CAP = self.get_CAP_from_file()

        numberOfRecordingSites = np.shape(CAP)[0]

        if numberOfRecordingSites > 2:

            numberOfPlots = min(maxNumberOfSubplots, numberOfRecordingSites-1)

            eletrodeSelection = np.floor(np.linspace(1,numberOfRecordingSites-1, numberOfPlots))

            # Subplots
            f, axarr = plt.subplots(numberOfPlots, sharex=True)

            for i in range(numberOfPlots):

                electrodeIndex = eletrodeSelection[i]

                CAPSingleElectrode =  CAP[electrodeIndex,:]
                distanceFromOrigin = self.saveParams['L']/numberOfRecordingSites*electrodeIndex

                axarr[i].plot(time, CAPSingleElectrode)
                axarr[i].set_title('distance ' + str(distanceFromOrigin) + ' [um]')
                axarr[i].set_ylabel('CAP [mV]')

                if i == numberOfPlots - 1:
                    axarr[i].set_xlabel('time [ms]')
        else:
            fig = plt.figure()
            CAPSingleElectrode =  CAP[numberOfRecordingSites-1,:]
            distanceFromOrigin = self.recording_elec_pos[0]

            plt.plot(time, CAPSingleElectrode)
            plt.title('distance ' + str(distanceFromOrigin) + ' [um]')
            plt.ylabel('CAP [mV]')
            plt.xlabel('time [ms]')

        plt.savefig(self.basePath+'CAP1D.png')


    def plot_CAP2D(self):

        # first load the desired data from file
        time, CAP = self.get_CAP_from_file()

        # check if plotting makes sense
        numberOfRecordingSites = np.shape(CAP)[0]

        if numberOfRecordingSites <= 10:
            print 'Plotting of the CAP in two dimensions (time, space) does not make sense with fewer than 10 electrodes. ' \
                  'Please select anotherplotting mechanism or restart the simulation with more electrodes.'
            return


        # print as an image
        fig = plt.figure()
        im = plt.imshow(CAP, cmap=plt.get_cmap('gist_stern'), interpolation='none', aspect='auto')#, norm=LogNorm(vmin=CAPmin, vmax=CAPmax))

        # correct xticks (from samples to ms)
        numberOfXTicks = 10
        tick_locs = np.round(np.linspace(0,np.shape(CAP)[1],numberOfXTicks))
        tick_lbls = np.round(np.linspace(0,self.saveParams['tStop'],numberOfXTicks))
        plt.xticks(tick_locs, tick_lbls, fontsize=12)

        # correct yticks (from electrodes to distance)
        numberOfYTicks = 10
        tick_locs = np.round(np.linspace(0,np.shape(CAP)[0],numberOfYTicks))
        tick_lbls = np.round(np.linspace(0,self.saveParams['L'],numberOfYTicks))
        plt.yticks(tick_locs, tick_lbls, fontsize=12)

        # add titles, axis labels and colorbar
        fig.suptitle('Compound action potential [uV] over space and time', fontsize=20)
        plt.xlabel('time [ms]')
        plt.ylabel('Distance from axon origin [um]')
        cbar = plt.colorbar(im)

        plt.savefig(self.basePath+'CAP2D.png')

    def loadVoltageFromFile(self):
        # get the whole CAP, can be signle electrode or multiple
        directory = getDirectoryName("V", self.basePath)
        try:
            newestFile = max(glob.iglob(directory+'*.[Dd][Aa][Tt]'), key=os.path.getctime)
        except ValueError:
            print 'No voltage calculation has been performed yet with this set of parameter.'
            return

        # load the raw voltage file
        timeStart = time.time()
        # Vraw = np.transpose(np.loadtxt(newestFile))
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

    def plot_voltage(self):

        timeRec, voltageMatrices = self.loadVoltageFromFile()

        # now plot
        numberOfAxons = np.shape(voltageMatrices)[0]
        numberOfPlots = min(6, numberOfAxons)

        axonSelection = np.floor(np.linspace(0,numberOfAxons-1, numberOfPlots))

        f, axarr = plt.subplots(numberOfPlots, sharex=True)

        # colors
        jet = plt.get_cmap('jet')

        # for colorbar first check what is the longest axon
        maxAxonLength = 0
        for i in range(len(axonSelection)):
            axonIndex = int(axonSelection[i])
            maxAxonLength = max(maxAxonLength, self.axons[axonIndex].L)

        cNorm = colors.Normalize(vmin=0, vmax=int(maxAxonLength)-1)#len(diameters_m)-1)#
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

        for i in range(len(axonSelection)):
            axonIndex = int(axonSelection[i])

            voltageMatrix = np.transpose(voltageMatrices[axonIndex])

            # find out whether axon is myelinated or not
            isMyelinated = (type(self.axons[axonIndex]) == Myelinated)

            axonDiameter = self.axons[axonIndex].fiberD
            currentNumberOfSegments = np.shape(voltageMatrix)[1]
            currentAxonLength = self.axons[axonIndex].L

            if not isMyelinated:
                for j in range(currentNumberOfSegments):
                    colorVal = scalarMap.to_rgba(int(j * currentAxonLength / currentNumberOfSegments))
                    axarr[i].plot(timeRec, voltageMatrix[:,j], color=colorVal)

                axarr[i].set_ylabel('Voltage [mV]')
                axarr[i].set_xlabel('time [ms]')
                axarr[i].set_title('Voltage of unmyelinated axon with diameter ' + str(axonDiameter) + ' um')
            else:
                Nnodes = self.axons[axonIndex].axonnodes

                numberOfRecordingSites = np.shape(voltageMatrix)[1]

                nodePositions = range(0,(Nnodes-1)*11,11)

                nodeDistance = self.axons[axonIndex].lengthOneCycle

                nodeCounter = 0
                for j in nodePositions:
                    colorVal = scalarMap.to_rgba(int(nodeCounter * nodeDistance))
                    axarr[i].plot(np.array(timeRec), np.array(voltageMatrix[:,j]), color=colorVal)
                    nodeCounter += 1

                axarr[i].set_ylabel('Voltage [mV]')
                # axarr[i].set_ylim([-100,100])
                axarr[i].set_xlabel('time [ms]')
                axarr[i].set_title('Voltage of nodes of myelinated axon with diameter ' + str(axonDiameter) + ' um')

            # make room for colorbar
            f.subplots_adjust(right=0.8)

            # add colorbar axis
            axColorbar = f.add_axes([0.85, 0.15, 0.05, 0.7])

            cb1 = mpl.colorbar.ColorbarBase(axColorbar, cmap=jet,
                                            norm=cNorm,
                                            orientation='vertical')
            cb1.set_label('length [um]')

        plt.savefig(self.basePath+'voltage.png')

    def create_axon(self, axonPosition):

        # first decide by chance whether myelinated or unmyelinated
        axonTypeIndex = np.random.choice(2,1,p = [self.p_A, self.p_C])
        axonTypes = ['m', 'u']
        axonType = axonTypes[axonTypeIndex]

        # then get diameter. Either drawn from distribution or constant.
        axonDiameter = self.getDiam(axonType)

        # axonCoords = np.row_stack((np.concatenate(([0], axonPosition)), np.concatenate(([self.bundleLength], axonPosition))))

        if True:
            # calculate the random axon coordinates
            axonCoords = createGeometry.create_random_axon(self.bundleCoords, self.radius_bundle, axonPosition,
                                                           self.segmentLengthAxon, randomDirectionComponent=self.randomDirectionComponent)
        else:
            axonCoords = np.column_stack(axonPosition, np.concatenate(([axonPosition[0] + self.bundleLength], axonPosition[1:2])))

        if axonTypeIndex == 1:
            unmyel = copy.copy(self.unmyelinated)
            unmyel['diam'] = axonDiameter
            axonParameters = dict( {'coord': axonCoords},**unmyel)
            #axonParameters = dict( {'coord': axonPosition},**unmyel)
            self.axons.append(Unmyelinated(**axonParameters))


        elif axonTypeIndex == 0:
            myel = copy.copy(self.myelinated_A)
            myel['fiberD'] = axonDiameter
            axonParameters = dict( {'coord':axonCoords},**myel)
            # axonParameters = dict( {'coord':axonPosition},**myel)
            self.axons.append(Myelinated(**axonParameters))
        else:
            "Error in the draw of the axon type!"

        self.axons[-1].axonPosition = axonPosition

        # self.stim = Stimulus(self.stim_type,self.axons[i], delay[i],self.stim_dur,self.amp, self.freq,self.duty_cycle, self.stim_coord, self.waveform)

    def simulateAxons(self):

        # where are the electrodes
        # [X,Y,Z,N] = self.setup_recording_elec()
        [X,Y,Z] = self.setup_recording_elec()

        for axonIndex in range(self.virtual_number_axons):

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

            ### LAUNCH SIMULATION ###
            temp = time.time()

            # create the neuron object specified in the axon class object
            axon.create_neuron_object()

            if not self.stim_type == 'NONE':
                # connect to stimulus
                self.stim.connectAxon(axon)

            try:
                # connect up stream nerve spiking
                self.upstreamSpiking.connectAxon(axon)
            except:
                pass

            if axonIndex == 0:
            # record time variable
                self.trec = h.Vector()
                self.trec.record(h._ref_t)



            axon.simulate()
            self.electrodes.append(refextelectrode.RecExtElectrode(axon, **electrodeParameters))
            elapsed1 = time.time()-temp
            print "Elapsed time calculate voltage and membrane current: " + str(elapsed1)

            temp = time.time()
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

            # keep axon objects clean
            for vrec in axon.vreclist:
                vrec = None
            axon.vreclist = None



    def store_geometry(self):
        self.geometry_parameters = [self.axons[0].xstart,self.axons[0].ystart,self.axons[0].zstart,self.axons[0].xend,self.axons[0].yend,self.axons[0].zend,self.axons[0].area,self.axons[0].diam,self.axons[0].length,self.axons[0].xmid,self.axons[0].ymid,self.axons[0].zmid]

    def save_electrode(self,i):
        directory = getDirectoryName("elec", self.basePath)

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

        directory = getDirectoryName("elec", self.basePath)
        filename = "electrode_"+str(elecIndex)+".dat"

        t0 = time.time()
        electrodeData = np.loadtxt(directory + filename, unpack=True)
        print "loaded electrode "+ str(elecIndex) +  " in " + str(time.time()-t0)

        return electrodeData

    def compute_CAP_fromfiles(self):
        temp = time.time()

        monopolar = len(self.recording_elec_pos) == 1

        # variable to save the sum over all axons
        self.sum_CAP = np.zeros((self.number_elecs,len(self.trec)))

        # variable to save the extracellular signal from each cell separately, at the last electrode position.
        self.AP_axonwise = np.zeros((self.virtual_number_axons, len(self.trec)))

        # load the recordings for every axon one by one and add them.
        for elecIndex in range(self.virtual_number_axons):
            electrodeData = self.load_one_electrode(elecIndex)

        # The contactpoints that constitute one cuff electrode ring have to be recovered, summed up together per
        # recording location along the axon
            for i in range(self.number_elecs):
                if monopolar:
                    contactPointIndices = range(self.number_contact_points*i, self.number_contact_points*(1+i))
                    sumOverContactPoints = np.sum(electrodeData[contactPointIndices, :], 0)
                else:
                    contactPointIndicesPole1 = range(self.number_contact_points*2*i, self.number_contact_points*(1+2*i))
                    contactPointIndicesPole2 = range(self.number_contact_points*(2*i+1), self.number_contact_points*(2*(i+1)))
                    sumOverContactPoints = np.sum(electrodeData[contactPointIndicesPole1, :] - electrodeData[contactPointIndicesPole2, :], 0)

                self.sum_CAP[i,:] = self.sum_CAP[i,:] +  sumOverContactPoints

                if i == self.number_elecs-1:
                    self.AP_axonwise[elecIndex,:] = sumOverContactPoints

        elapsed = time.time()-temp
        print "Elapsed time to compute CAP " + str(elapsed) + " \n"


    def setup_recording_elec(self):

        if True:
            # calculte recording electrode positions for a 3D shaped bundle
            electrodePositions = createGeometry.electrodePositionsBundleGuided(self.bundleCoords, self.radius_bundle,
                                                                               self.number_elecs, self.number_contact_points,
                                                                               self.recording_elec_pos)
            X, Y, Z = electrodePositions[:,0], electrodePositions[:,1], electrodePositions[:,2]
        else:

            if (self.number_elecs == 1):
                # if one recording site
                if len(self.recording_elec_pos) == 1:
                    # if monopolar
                    X = np.zeros(self.number_contact_points)+self.recording_elec_pos
                elif len(self.recording_elec_pos) == 2:
                    # if bipolar
                    X = np.repeat(self.recording_elec_pos,self.number_contact_points,axis=0)
                elif len(self.recording_elec_pos) > 2 or len(self.recording_elec_pos) == 0:
                    # if wrong number of poles entered
                    raise TypeError("Only monopolar and bipolar recording are supported")
            else:
                if len(self.recording_elec_pos) >1:
                    raise TypeError("Please use only the monopolar configuration of 'recording_elec' to record along the bundle")
                X1 = [np.linspace(0, self.recording_elec_pos[0], self.number_elecs)]
                X = np.repeat(X1,self.number_contact_points, axis=0)
                X = X.flatten()
            angles = np.linspace(0,360, self.number_contact_points, endpoint = False)
            Y1 = np.round(self.radius_bundle*np.cos(angles*np.pi/180.0),2)
            Z1 = np.round(self.radius_bundle*np.sin(angles*np.pi/180.0),2)
            Y = np.tile(Y1,self.number_elecs*len(self.recording_elec_pos))
            Z = np.tile(Z1,self.number_elecs*len(self.recording_elec_pos))
            N = np.empty((self.number_contact_points*self.number_elecs*len(self.recording_elec_pos), 3))
            for i in xrange(N.shape[0]):
                N[i,] = [1, 0, 0] #normal vec. of contacts

        self.electrodeCoords = np.transpose(np.row_stack((X,Y,Z)))

        return [X,Y,Z]#,N]

    def build_disk(self,number_of_axons,radius_bundle):
        ### AXONS POSITIONS ON THE DISK ###
        """
        Partly from http://blog.marmakoide.org/?p=1
        """
        n = number_of_axons

        radius = radius_bundle* np.sqrt(np.arange(n) / float(n))

        golden_angle = np.pi * (3 - np.sqrt(5))
        theta = golden_angle * np.arange(n)

        self.axons_pos = np.zeros((n, 2))
        self.axons_pos[:,0] = np.cos(theta)
        self.axons_pos[:,1] = np.sin(theta)
        self.axons_pos *= radius.reshape((n, 1))

    def get_filename(self, recordingType):

        directory = getDirectoryName(recordingType, self.basePath)
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

    def get_filename_for_draw(self):
        self.filename = "p_A"+str(self.p_A)+"_p_C"+str(self.p_C)+'nb_axons'+str(self.number_of_axons)+'bundle_radius'+str(self.radius_bundle)+'.dat'

        return self.filename

    def getDiam(self, axonType):

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
            givenDiameter = self.unmyelinated['diam']

            if isinstance(givenDiameter, float) or isinstance(givenDiameter, int):
                return givenDiameter
            elif isinstance(givenDiameter, dict):
                fiberD = self.unmyelinated['diam']
                densities = fiberD['densities']

                sum_d = float(sum(densities))
                normalize_densities = [x / sum_d for x in densities]

                draw_diam = np.random.choice(len(normalize_densities),1,p = normalize_densities)
                D = fiberD['diameters'][draw_diam]
                diam = D
            else:
                raise('Something is wrong with your given axon diameter for unmyelinated axons.')


        else:
            raise('Wrong axon type given to function getDiam. Valid ones: u or m')

        return diam
