from classes import *
from classes_plot import *
import time
time.sleep(0) # delays for x seconds

""" McIntyre Model parameters
celsius=37			
v_init=-80 //mV//  		
dt=0.005 //ms//         	
tstop=10"""
h.celsius = 33 # set temperature in celsius
h.tstop = 30 # set simulation duration (ms)
h.dt = 0.0025 # set time step (ms)
h.finitialize(-65) # initialize voltage state

# Set parameters
calculationFlag = False
plottingFlag = True

plotCAP = False
plotCAP1D = False
plotCAP2D = False

plotVoltage = True

# bundle characteristics
p_A = [1.]#[0.175,0.1,1.0, 0.0] # share of myelinated fibers
fiberD_A = 16.0 #um diameter myelinated axons 'draw' OR one of 5.7, 7.3, 8.7, 10.0, 11.5, 12.8, 14.0, 15.0, 16.0
fiberD_C = 1.5#'draw'


radius_bundle = 150.0 #um Radius of the bundle (typically 0.5-1.5mm)
draw_distribution = True #Boolean stating the distribution of fibre should be drawn
number_of_axons =  50
lengthOfBundle = 1000


# stimulus characteristics
stim_types = ["INTRA"]#, "INTRA", "EXTRA"]
waveforms = ["MONOPHASIC"]#,"MONOPHASIC", "BIPHASIC"]
frequencies = [0.1]#,0.1,0.1]
duty_cycles = [0.01]#[0.001]#,0.01,0.005]
amplitudes = [4.0]#,2.0,0.5]
stimDur = [10]

# recoding params
number_contact_points=  8 #Number of points on the circle constituing the cuff electrode
recording_elec_pos = [1000] #[10000], #Position of the recording electrode along axon in um, in "BIPOLAR" case the position along axons should be given as a couple [x1,x2]
number_elecs =  150#150, #number of electrodes along the bundle

# Do not change from here


### !!! USING MYELINATED DISTRUTION NOT COMPATIBLE WITH VOLTAGE RECORDING ###
### ANYWAY VOLTAGE RECORDING SHOULD LOGICALLY BE PERFORMED FOR ONLY ONE AXON AT A TIME ###
myelinatedDistribution = {
    'densities':[100,300,1150,2750,3650,2850,1750,900,500,250,200,150,110,100,110,100,105], #fibers densities can be given either in No/mm2 or percentage
    'diameters': [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7.,7.5, 8., 8.5, 9.],  # corresponding diameters for each densities
}

unmyelinatedDistribution = {
    'densities':[250,1250,5000,8000,9800,10200,8900,7600,5700,4000,3900,2300,2000,1300,900,750,600,600,500,250], #fibers densities can be given either in No/mm2 or percentage
    'diameters': [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.], # corresponding diameters for each densities
}

if fiberD_A == 'draw':
    del fiberD_A
    fiberD_A = myelinatedDistribution
# else:
#     ### PARAMETERS BELOW ARE ACTUALLY SET IN THE CODE BUT USEFUL FOR THE HEADER AND PLOTS ###
#     ## Values of nodelength and paralength1 are constant
#     ## Values of paralength2 and interlength are to be set according to the chosen fiberD value
#     nodelength = 1.0 #um
#     paralength1 = 1.3 #um
#     [paralength2_A, interlength_A] = fiberD_dependent_param(fiberD_A, nodelength, paralength1)
if fiberD_C == 'draw':
    del fiberD_C
    fiberD_C = unmyelinatedDistribution



for VoltCAPSelector in [2]:#[1,2]:
    rec_CAP = (VoltCAPSelector==1)
    rec_v = (VoltCAPSelector==2)
    for j in range(len(duty_cycles)):
        for k in range(len(p_A)):#[0]:#range(1,len(p_A)):#
            stimulusParameters = {
                'jitter_para': [0,0], #Mean and standard deviation of the delay
                'stim_type': stim_types[j], #Stimulation type either "INTRA" or "EXTRA"
                # stim_coord is NOT USED default value is used directly
                'stim_coord': [[0,50,0]], # spatial coordinates  of the stimulating electrodes, example for tripolar case=[[xe0,ye0,ze0], [xe1,ye1,ze1], [xe2,ye2,ze2]] (order is important with the middle being the cathode), INTRA case only use the position along x for IClamp
                'amplitude': amplitudes[j], # Pulse amplitude (nA)
                'freq': frequencies[j], # Frequency of the sin pulse (kHz)
                'duty_cycle': duty_cycles[j], # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
                'stim_dur' : stimDur[j], # Stimulus duration (ms)
                'waveform': waveforms[j], # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
    }

            recordingParameters = {
                "number_contact_points": number_contact_points, #Number of points on the circle constituing the cuff electrode
                'recording_elec_pos': recording_elec_pos,#[10000], #Position of the recording electrode along axon in um, in "BIPOLAR" case the position along axons should be given as a couple [x1,x2]
                'number_elecs': number_elecs,#150, #number of electrodes along the bundle
                'dur': h.tstop, # Simulation duration (ms)
                'rec_CAP': rec_CAP, #If false means we avoid spending time using LFPy functions
            }
            myelinatedParametersA = {
                'name': "myelinated_axonA", # axon name (for neuron)
                'Nnodes': 11, #Number of nodes
                'fiberD': fiberD_A, #fiberD_A, #Diameter of the fiber
                'layout3D': "DEFINE_SHAPE", # either "DEFINE_SHAPE" or "PT3D" using hoc corresponding function
                'rec_v': rec_v, # set voltage recorders True or False
                # 'nodelength' : nodelength, #set node length (um)
                # 'paralength1': paralength1, #set the length of the nearest paranode segment to the node
                # 'paralength2': paralength2_A,  #set the length of the second paranode segment followed by the internodes segments
                # 'interlength': interlength_A, #set the length of the internode part comprising the 6 segments between two paranodes2
            }


            unmyelinatedParameters = {
                'name': "unmyelinated_axon", # axon name (for neuron)
                'L': lengthOfBundle,#Axon length (micrometer)
                'diam': fiberD_C, #Axon diameter (micrometer)
                'cm' : 1.0, #Specific membrane capacitance (microfarad/cm2)
                'Ra': 200.0, #Specific axial resistance (Ohm cm)
                'layout3D': "DEFINE_SHAPE", # either "DEFINE_SHAPE" or "PT3D" using hoc corresponding function
                'rec_v': rec_v, # set voltage recorders True or False
            }
            bundleParameters = {         #parameters for Bundle classe
                'radius_bundle': radius_bundle, #um Radius of the bundle (typically 0.5-1.5mm)
                'draw_distribution': draw_distribution, #Boolean stating the distribution of fibre should be drawn
                'number_of_axons': number_of_axons,#640, # Number of axons in the bundle
                'p_A': p_A[k], # Percentage of myelinated fiber type A
                'p_C': 1-p_A[k], #Percentage of unmyelinated fiber type C
                'myelinated_A': myelinatedParametersA, #parameters for fiber type A
                'unmyelinated': unmyelinatedParameters, #parameters for fiber type C
            }


            Parameters1 = dict(bundleParameters, **stimulusParameters)
            Parameters = dict(Parameters1, **recordingParameters)

            if calculationFlag:

                bundle = Bundle(**Parameters)

                if rec_CAP:
                    save_CAP_tofile(bundle,Parameters)
                if rec_v:
                    # When saving voltage to file limit the number of axons to 10. If 100 unmyelinated it produces a 1Go file, and if 100 myelinated 2Go.
                    save_voltage_tofile(bundle,Parameters)

                bundle = None

            if plottingFlag:


                saveParams={'elecCount': len(recording_elec_pos), 'dt': h.dt, 'tStop': h.tstop, 'p_A': bundleParameters['p_A'],
                    'myelinatedDiam': myelinatedParametersA['fiberD'], 'unmyelinatedDiam': unmyelinatedParameters['diam'],
                    'L': unmyelinatedParameters['L'], 'stimType': stimulusParameters['stim_type'], 'stimWaveform' : stimulusParameters['waveform'],
                    'stimDutyCycle': stimulusParameters['duty_cycle'], 'stimAmplitude' : stimulusParameters['amplitude']}

                if plotCAP:

                    # get the whole CAP, can be signle electrode or multiple
                    directory = getDirectoryName("CAP", **saveParams)
                    try:
                        newestFile = max(glob.iglob(directory+'*.[Dd][Aa][Tt]'), key=os.path.getctime)
                    except ValueError:
                        print 'No CAP calculation has been performed yet with this set of parameter.'
                        quit()

                    CAPraw = np.transpose(np.loadtxt(newestFile))
                    time = CAPraw[0,:]
                    CAP = CAPraw[1:,:]

                    if plotCAP1D:

                        numberOfRecordingSites = np.shape(CAP)[0]
                        numberOfPlots = min(5, numberOfRecordingSites)

                        eletrodeSelection = np.floor(np.linspace(0,numberOfRecordingSites-1, numberOfPlots))

                        # Subplots
                        f, axarr = plt.subplots(numberOfPlots, sharex=True)

                        for i in range(numberOfPlots):

                            electrodeIndex = eletrodeSelection[i]

                            CAPSingleElectrode =  CAP[electrodeIndex,:]
                            distanceFromOrigin = saveParams['L']/numberOfRecordingSites*electrodeIndex

                            axarr[i].plot(CAPSingleElectrode)
                            axarr[i].set_title('distance ' + str(distanceFromOrigin) + ' [um]')
                            axarr[i].set_ylabel('CAP [uV]')

                            if i == numberOfPlots - 1:
                                axarr[i].set_xlabel('time [ms]')




                    if plotCAP2D:

                        # print as an image
                        fig = pyl.figure()
                        im = pyl.imshow(CAP, cmap=pyl.get_cmap('gist_stern'), interpolation='none', aspect='auto')#, norm=LogNorm(vmin=CAPmin, vmax=CAPmax))

                        # correct xticks (from samples to ms)
                        numberOfXTicks = 10
                        tick_locs = np.round(np.linspace(0,np.shape(CAP)[1],numberOfXTicks))
                        tick_lbls = np.round(np.linspace(0,saveParams['tStop'],numberOfXTicks))
                        plt.xticks(tick_locs, tick_lbls, fontsize=12)

                        # correct yticks (from electrodes to distance)
                        numberOfYTicks = 10
                        tick_locs = np.round(np.linspace(0,np.shape(CAP)[0],numberOfYTicks))
                        tick_lbls = np.round(np.linspace(0,saveParams['L'],numberOfYTicks))
                        plt.yticks(tick_locs, tick_lbls, fontsize=12)

                        # add titles, axis labels and colorbar
                        fig.suptitle('Compound action potential [uV] over space and time', fontsize=20)
                        pyl.xlabel('time [ms]')
                        pyl.ylabel('Distance from axon origin [um]')
                        cbar = pyl.colorbar(im)

                    # finally show result
                    pyl.show()


                if plotVoltage:

                    # get the whole CAP, can be signle electrode or multiple
                    directory = getDirectoryName("V", **saveParams)
                    try:
                        newestFile = max(glob.iglob(directory+'*.[Dd][Aa][Tt]'), key=os.path.getctime)
                    except ValueError:
                        print 'No voltage calculation has been performed yet with this set of parameter.'
                        quit()

                    # load the raw voltage file
                    Vraw = np.transpose(np.loadtxt(newestFile))

                    time = Vraw[0,1:] # extract time vector
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
                    numberOfPlots = 2#min(5, numberOfAxons)

                    axonSelection = np.floor(np.linspace(0,numberOfAxons-1, numberOfPlots))

                    f, axarr = plt.subplots(numberOfPlots, sharex=True)

                    # colors
                    jet = plt.get_cmap('jet')

                    for i in range(len(axonSelection)):
                        voltageMatrix = np.transpose(voltageMatrices[int(axonSelection[i])])

                        currentNumberOfSegments = np.shape(voltageMatrix)[1]

                        cNorm = colors.Normalize(vmin=0, vmax=currentNumberOfSegments-1)#len(diameters_m)-1)#
                        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

                        for j in range(currentNumberOfSegments):
                            colorVal = scalarMap.to_rgba(j)
                            axarr[i].plot(voltageMatrix[:,j], color=colorVal)
                        axarr[i].set_title('distance ' + str(555) + ' [um]')
                        axarr[i].set_ylabel('Voltage [mV]')

                    # finally show result
                    pyl.show()
                    print 'plotit'




