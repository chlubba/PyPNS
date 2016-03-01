from axonClass import *
from upstreamSpikingClass import *
from stimulusClass import *

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
from mpl_toolkits.mplot3d import Axes3D



from nameSetters import getDirectoryName, getFileName

class Bundle(object):
    """
    radius_bundle: Radius of the bundle (typically 0.5-1.5mm)
    number_of_axons: Number of axons in the bundle
    p_A: Percentage of myelinated fiber type A
    p_B: Percentage of myelinated fiber type B
    p_C: Percentage of unmyelinated fiber type C
    number_contact_points: Number of points on the circle constituing the cuff electrode
    recording_elec: Position of the recording electrode along axon in um, in "BIPOLAR" case the position along axons should be given as a couple [x1,x2]
    jitter_para: Mean and standard deviation of the jitter
    stim_type: Stimulation type either "INTRA" or "EXTRA" for INTRA/EXTRA_cellular stimulation
    stim_coord: spatial coordinates  of the stimulating electrode, example for tripolar case=[[xe0,ye0,ze0], [xe1,ye1,ze1], [xe2,ye2,ze2]] (order is important with the middle being the cathode)
    amplitude: pulse amplitude (nA)
    freq: frequency of the sin pulse (kHz)
    duty_cycle: Percentage stimulus is ON for one cycl
    stim_dur : stimulus duration (ms)
    dur: simulation duration (ms)
    waveform: Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
    number_of_elecs: number of electrodes along the bundle
    recording_type: "BIPOLAR" or "MONOPOLAR"
    myelinated_A: parameters for fiber type A
    umyelinated:  parameters for fiber type C
    draw_distribution: Boolean stating the distribution of fibre should be drawn


    """
    def __init__(self, radius_bundle, draw_distribution, number_of_axons, p_A, p_C, number_contact_points, recording_elec_pos, jitter_para, stim_type, stim_coord, duty_cycle, freq, amplitude, stim_dur, dur, number_elecs, myelinated_A, unmyelinated,rec_CAP, waveform):

        self.draw_distribution = draw_distribution
        self.myelinated_A =  myelinated_A
        self.unmyelinated =  unmyelinated

        self.waveform = waveform
        self.stim_type = stim_type
        #self.stim_coord = stim_coord#[0,500,0] # [xe,ye,ze] #um
        self.stim_coord = [[0, radius_bundle*math.cos(math.pi/number_contact_points), radius_bundle*math.sin(math.pi/number_contact_points)]]
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
        self.radius_bundle = radius_bundle # um
        self.electrodes = []
        self.voltages = []
        self.number_contact_points = number_contact_points
        self.number_elecs = number_elecs
        self.recording_elec_pos = recording_elec_pos #um

        self.build_disk(self.number_of_axons,self.radius_bundle)

        self.saveParams={'elecCount': len(self.recording_elec_pos), 'dt': h.dt, 'tStop': h.tstop, 'p_A': self.p_A,
                    'myelinatedDiam': self.myelinated_A['fiberD'], 'unmyelinatedDiam': self.unmyelinated['diam'],
                    'L': self.unmyelinated['L'], 'stimType': self.stim_type, 'stimWaveform' : self.waveform,
                    'stimDutyCycle': self.duty_cycle, 'stimAmplitude' : self.amp}

        # ### JITTER (random gaussian delay for individual fiber stimulation) ###
        # self.delay_mu, self.delay_sigma = jitter_para[0], jitter_para[1] # mean and standard deviation
        # if (jitter_para[0] == 0) & (jitter_para[1] == 0):
        #     delay = np.zeros(self.number_of_axons)
        # else:
        #     delay = np.random.normal(self.delay_mu, self.delay_sigma, self.number_of_axons)

        # create axons within wanted nerve-slice


        # create axons within one fourth slice of the whole bundle.
        self.virtual_number_axons = 0
        for i in range(self.number_of_axons):
            # if within the desired subsection of the bundle, create axon
            if True:#((self.axons_pos[i,1]>=0 and self.axons_pos[i,1]< self.axons_pos[i,0]) or (self.axons_pos[i,0]== 0 and self.axons_pos[i,1]==0)):
                # print self.axons_pos[i,:]
                self.create_axon(self.axons_pos[i,:])
                self.virtual_number_axons +=1
                print "Number axons created:" + str(self.virtual_number_axons)


        if not self.stim_type == 'NONE':
            # create Simulus instace used for all axons
            self.stim = Stimulus(self.stim_type, self.stim_dur,self.amp, self.freq,self.duty_cycle, self.stim_coord, self.waveform)


        # # create upstream activity
        # self.upstreamSpiking = UpstreamSpiking(self.number_of_axons, tStart=0., tStop=h.tstop, lambd = 1000.)


        # # connect to stimulus
        # for i in range(self.virtual_number_axons):
        #     axon = self.axons[i]
        #     self.stim.connectAxon(axon)

    def addUpstreamSpiking(self, tStart=0., tStop=h.tstop, lambd = 1000., correlation = 0.1):
        # create upstream activity
        self.upstreamSpiking = UpstreamSpiking(self.number_of_axons, tStart=tStart, tStop=tStop, lambd=lambd, correlation=correlation)

    def simulateBundle(self):

        self.simulateAxons()

        if self.rec_CAP:
            self.compute_CAP_fromfiles()
            self.save_CAP_to_file()
            self.sum_CAP = None

        self.save_voltage_to_file()

        # get rid of the all Neuron objects to be able to pickle the bundle-class. pickle. lick it.
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
            except AttributeError:
                pass
            try:
                for vrec in axon.vreclist:
                    vrec = None
                axon.vreclist = None
            except AttributeError:
                pass
            axon.allseclist = None
            # also delete unnecessary data that will no longer been used to keep the pickled file small
            axon.imem = None
        self.voltages = None



    def save_CAP_to_file(self):
        # filename = self.get_filename("CAP")
        filename = getFileName("CAP", self.saveParams)
        print "Save location for CAP file: " + filename

        # maybe add the header later. Right now we assume that the simulation is defined by the bundle object that get
        # always generated during the whole simulation. If files are opened independently from a bundle object, such a
        # header would be useful.
        # header = repr(parameters)
        DataOut = np.array(self.trec)

        for i in range(len(self.sum_CAP)):
                DataOut = np.column_stack( (DataOut, np.array(self.sum_CAP[i])))

        np.savetxt(filename, DataOut)

    def save_voltage_to_file(self):
        # filename = self.get_filename("V")
        filename = getFileName("V", self.saveParams)
        print "Save location for voltage file: " + filename
        #header= repr(parameters)

        DataOut = np.concatenate(([0],np.array(self.trec)))

        voltages = np.array(self.voltages)

        if np.size(voltages) == 0:
            return

        # as np.array(...) only converts outer Vector to python-readable format, we need to iterate through elements to convert
        # inner vectors where the actual voltage signals are stored.

        for i in range(len(voltages)):
            voltageSingleAxon = np.transpose(np.array(voltages[i]))

            # append the sectionlength in the first column in order to differentiate different axons later
            numberOfSegments = np.shape(voltageSingleAxon)[1]
            numberOfSegmentsArray = np.multiply(np.ones(numberOfSegments),np.array(numberOfSegments))
            voltageSingleAxonFormatted = np.row_stack((numberOfSegmentsArray, voltageSingleAxon))

            DataOut = np.column_stack( (DataOut, voltageSingleAxonFormatted))
        np.savetxt(filename, DataOut)#, header=header)

    def get_CAP_from_file(self):

        # get the whole CAP, can be single electrode or multiple
        directory = getDirectoryName("CAP", **self.saveParams)
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
            print 'Bundle has not been run yet. No geometry information was generated in Neuron.'
            return

        fig = plt.figure()
        ax = fig.gca(projection='3d')

        axonID = 0
        for axon in self.axons:
            if type(axon) == Myelinated:
                style = '-'
            else:
                style = '--'
            ax.plot(axon.xmid, axon.ymid, axon.zmid, style, label='axon '+str(axonID))
            axonID += 1
        plt.legend()

        elecCoords = self.electrodeCoords
        ax.scatter(elecCoords[:,0], elecCoords[:,1], elecCoords[:,2])
        for i in range(self.number_elecs):
            selectionIndices = range(i, self.number_contact_points*self.number_elecs, self.number_elecs)
            ringCoords = elecCoords[selectionIndices,:]
            ringCoords = np.row_stack((ringCoords, ringCoords[0,:]))
            ax.plot(ringCoords[:,0], ringCoords[:,1], ringCoords[:,2], color=[0.8,0.8,0.8])

    def plot_CAP1D_singleAxon(self, axonID):

        CAP = self.CAP

        print 'hm'


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

    def plot_voltage(self):
        # get the whole CAP, can be signle electrode or multiple
        directory = getDirectoryName("V", **self.saveParams)
        try:
            newestFile = max(glob.iglob(directory+'*.[Dd][Aa][Tt]'), key=os.path.getctime)
        except ValueError:
            print 'No voltage calculation has been performed yet with this set of parameter.'
            quit()

        # load the raw voltage file
        timeStart = time.time()
        Vraw = np.transpose(np.loadtxt(newestFile))
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

        # now plot
        numberOfAxons = np.shape(voltageMatrices)[0]
        numberOfPlots = min(6, numberOfAxons)

        axonSelection = np.floor(np.linspace(0,numberOfAxons-1, numberOfPlots))

        f, axarr = plt.subplots(numberOfPlots, sharex=True)

        # colors
        jet = plt.get_cmap('jet')

        for i in range(len(axonSelection)):
            axonIndex = int(axonSelection[i])

            voltageMatrix = np.transpose(voltageMatrices[axonIndex])

            # find out whether axon is myelinated or not
            isMyelinated = (type(self.axons[axonIndex]) == Myelinated)

            axonDiameter = self.axons[axonIndex].fiberD
            currentNumberOfSegments = np.shape(voltageMatrix)[1]

            if not isMyelinated:
                cNorm = colors.Normalize(vmin=0, vmax=currentNumberOfSegments-1)#len(diameters_m)-1)#
                scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)
                for j in range(currentNumberOfSegments):
                    colorVal = scalarMap.to_rgba(j)
                    axarr[i].plot(timeRec, voltageMatrix[:,j], color=colorVal)
                # axarr[i].set_title('distance ' + str(555) + ' [um]')
                axarr[i].set_ylabel('Voltage [mV]')
                axarr[i].set_xlabel('time [ms]')
                axarr[i].set_title('Voltage of unmyelinated axon with diameter ' + str(axonDiameter) + ' um')
            else:
                Nnodes = self.myelinated_A['Nnodes']

                numberOfRecordingSites = np.shape(voltageMatrix)[1]

                cNorm = colors.Normalize(vmin=0, vmax=numberOfRecordingSites-1)#len(diameters_m)-1)#
                scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

                nodePositions = range(0,(Nnodes-1)*11,11)

                # numberOfRecordingSites = np.shape(voltageMatrix)[1]
                #
                # cNorm = colors.Normalize(vmin=0, vmax=numberOfRecordingSites)#len(diameters_m)-1)#
                # scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)
                #
                # nodePositions = range(numberOfRecordingSites)

                for j in nodePositions:
                    colorVal = scalarMap.to_rgba(j)
                    axarr[i].plot(np.array(timeRec), np.array(voltageMatrix[:,j]), color=colorVal)

                axarr[i].set_ylabel('Voltage [mV]')
                axarr[i].set_xlabel('time [ms]')
                axarr[i].set_title('Voltage of nodes of myelinated axon with diameter ' + str(axonDiameter) + ' um')

    def create_axon(self, axonPosition):

        # first decide by chance whether myelinated or unmyelinated
        axonTypeIndex = np.random.choice(2,1,p = [self.p_A, self.p_C])
        axonTypes = ['m', 'u']
        axonType = axonTypes[axonTypeIndex]

        # then get diameter. Either drawn from distribution or constant.
        axonDiameter = self.getDiam(axonType)

        if axonTypeIndex == 1:
            unmyel = copy.copy(self.unmyelinated)
            unmyel['diam'] = axonDiameter
            axonParameters = dict( {'coord': axonPosition},**unmyel)
            self.axons.append(Unmyelinated(**axonParameters))


        elif axonTypeIndex == 0:
            myel = copy.copy(self.myelinated_A)
            myel['fiberD'] = axonDiameter
            axonParameters = dict( {'coord':axonPosition},**myel)
            self.axons.append(Myelinated(**axonParameters))
        else:
            "Error in the draw of the axon type!"

        # self.stim = Stimulus(self.stim_type,self.axons[i], delay[i],self.stim_dur,self.amp, self.freq,self.duty_cycle, self.stim_coord, self.waveform)

    def simulateAxons(self):

        # where are the electrodes
        [angles,X,Y,Z,N] = self.setup_recording_elec()

        runFlag = False
        for axonIndex in range(self.virtual_number_axons):

            axon = self.axons[axonIndex]

            # where is the axon
            axonPosition = axon.coord

            electrodeParameters = {         #parameters for RecExtElectrode class
                    'sigma' : 0.3,              #Extracellular potential
                    'x' : X,  #Coordinates of electrode contacts
                    'y' : Y-axonPosition[0],
                    'z' : Z-axonPosition[1],
                    'n' : 20,
                    'r' : 10,
                    'N' : N,
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


            if self.rec_CAP:
                axon.simulate()
                self.electrodes.append(refextelectrode.RecExtElectrode(axon, **electrodeParameters))
                elapsed1 = time.time()-temp
                print "Elapsed time to simulate CAP: " + str(elapsed1)
                temp = time.time()
                self.electrodes[axonIndex].calc_lfp()
                elapsed2 = time.time()-temp
                print "Elapsed time to calc lfp:" + str(elapsed2)
                self.save_electrode(axonIndex)
                self.electrodes[axonIndex]= None
                self.CAP_to_file = True

                # test if voltages can be recorded on the side
                self.voltages.append(axon.vreclist)

            elif axon.rec_v:
                # just run a normal NEURON simulation to record the voltage
                axon.set_voltage_recorders()
                h.run() #runFlag = True #
                self.voltages.append(axon.vreclist)
                elapsed = time.time()-temp
                print "Elapsed time to simulate V: " + str(elapsed)
            else:
                print "You are probably recording nothing for this axon"
                h.run() #runFlag = True #

            # delete the object
            axon.delete_neuron_object()

        # if runFlag:
        #     startTime = time.time()
        #     h.run()
        #     print 'Elapsed time for voltage simulation ' + str(time.time() - startTime)
            # ### DELETE THE PYTHON REFERENCE TO THE AXON TO DELETE THE AXON ###
            # for sec in h.allsec():
            #     for seg in sec:
            #         seg = None
            #     sec = None
            # if draw[i] == 1:
            #     self.axons[i].axon = None
            # else:
            #     self.axons[i].nodes = None
            #     self.axons[i].FLUTs = None
            #     self.axons[i].MYSAs = None
            #     self.axons[i].STINs = None


    def store_geometry(self):
        self.geometry_parameters = [self.axons[0].xstart,self.axons[0].ystart,self.axons[0].zstart,self.axons[0].xend,self.axons[0].yend,self.axons[0].zend,self.axons[0].area,self.axons[0].diam,self.axons[0].length,self.axons[0].xmid,self.axons[0].ymid,self.axons[0].zmid]

    def save_electrode(self,i):
        directory = getDirectoryName("elec", **self.saveParams)

        print "saving electrode: "+str(i)
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

        directory = getDirectoryName("elec", **self.saveParams)
        filename = "electrode_"+str(elecIndex)+".dat"

        t0 = time.time()
        electrodeData = np.loadtxt(directory + filename, unpack=True)
        print "loaded electrode "+ str(elecIndex) +  " in " + str(time.time()-t0)

        return electrodeData

    def compute_CAP_fromfiles(self):
        # directory = getDirectoryName("elec", **self.saveParams)
        #
        # temp = time.time()
        # CAP = []
        # print "loading electrode"
        # t0 = time.time()
        #
        # filename = "electrode_"+str(0)+".dat"
        # electrodesData = np.loadtxt(directory + filename, unpack=True)

        temp = time.time()

        self.sum_CAP = np.zeros((self.number_elecs,len(self.trec)))

        # load the recordings for every axon one by one and add them.
        for elecIndex in range(self.virtual_number_axons):
            electrodeData = self.load_one_electrode(elecIndex)

        # The contactpoints that constitute one cuff electrode ring have to be recovered, summed up together per
        # recording location along the axon
        # for i in range(self.number_contact_points):
        #     for j in range(self.number_elecs):
        #         CAP.append(electrodesData[i*self.number_elecs+j])

            for i in range(self.number_elecs):
                contactPointIndices = range(i*self.number_contact_points, (i+1)*self.number_contact_points)
                self.sum_CAP[i,:] = self.sum_CAP[i,:] +  np.sum(electrodeData[contactPointIndices, :], 0)

        # for k in range(1,self.virtual_number_axons):
        #     t0 = time.time()
        #     filename = "electrode_"+str(k)+".dat"
        #     electrodesData = np.loadtxt(directory + filename, unpack=True)
        #     print "electrode_"+str(k)+ "loaded in: "+str(time.time()-t0)
        #     for j in range(self.number_elecs):
        #         for i in range(self.number_contact_points):
        #             CAP[i*self.number_elecs+j] += electrodesData[i*self.number_elecs+j]
        #             self.sum_CAP[j,:] += CAP[i*self.number_elecs+j]
        #     del electrodesData

        elapsed = time.time()-temp
        print "Elapsed time to compute CAP:" + str(elapsed)


    def setup_recording_elec(self):
        if (self.number_elecs == 1):
            if len(self.recording_elec_pos) == 1:
                X = np.zeros(self.number_contact_points)+self.recording_elec_pos
            elif len(self.recording_elec_pos) == 2:
                X = np.repeat(self.recording_elec_pos,self.number_contact_points,axis=0)
            elif len(self.recording_elec_pos) > 2 or len(self.recording_elec_pos) == 0:
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
        Y = np.repeat(Y1,self.number_elecs*len(self.recording_elec_pos))
        Z = np.repeat(Z1,self.number_elecs*len(self.recording_elec_pos))
        N = np.empty((self.number_contact_points*self.number_elecs*len(self.recording_elec_pos), 3))
        for i in xrange(N.shape[0]):
            N[i,] = [1, 0, 0] #normal vec. of contacts

        self.electrodeCoords = np.transpose(np.row_stack((X,Y,Z)))

        return [angles,X,Y,Z,N]

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

        directory = getDirectoryName(recordingType, **self.saveParams)
        if not os.path.exists(directory):
            os.makedirs(directory)

        # filename = 'recording.dat'
        filename = recordingType+'.dat'

        number = 0
        filenameTemp = filename
        while os.path.isfile(directory+filenameTemp):
            number += 1
            print "Be careful this file name already exist ! We concatenate a number to the name to avoid erasing your previous file."
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
