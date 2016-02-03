from classes import *
from classes_plot import *


""" McIntyre Model parameters
celsius=37			
v_init=-80 //mV//  		
dt=0.005 //ms//         	
tstop=10"""
h.celsius = 33 # set temperature in celsius
h.tstop = 3e1 # set simulation duration (ms)
h.dt = 0.0025 # set time step (ms)
h.finitialize(-65) # initialize voltage state



#fiberD from 5.7, 7.3, 8.7, 10.0, 11.5, 12.8, 14.0, 15.0, 16.0
fiberD_A = 11.5 #um

### PARAMETERS BELOW ARE ACTUALLY SET IN THE CODE BUT USEFUL FOR THE HEADER AND PLOTS ###
## Values of nodelength and paralength1 are constant
## Values of paralength2 and interlength are to be set according to the chosen fiberD value
nodelength = 1.0 #um
paralength1 = 1.3 #um
[paralength2, interlength] = fiberD_dependent_param(fiberD_A, nodelength, paralength1)

Nnodes = 41
index_node = Nnodes-31
p_A = [1]
pos1 = nodelength*(Nnodes-2-index_node)+nodelength*0.5+2*(Nnodes-2-index_node)*(paralength1+paralength2)+6*(Nnodes-2-index_node)*interlength
pos2= pos1+paralength1/2.0
pos3 = pos2 + paralength1/2.0 + paralength2/2.0
pos4 = pos3 + paralength2/2.0 + 3*interlength
recording_pos = [pos1, pos2, pos3, pos4]
frequency = [0.1]
duty = [0.01]
stim_type = ["INTRA"]
amplitudes = [0.0]#[2.0]
waveform = ["MONOPHASIC"]
for i in range(len(duty)):
    for w in range(4):
        for k in range(1):
            for j in range(1):
                stimulusParameters = {
                    'jitter_para': [0,0], #Mean and standard deviation of the delay
                    'stim_type': stim_type[i], #Stimulation type either "INTRA" or "EXTRA" 
                    'stim_coord': [[0,50,0]], # spatial coordinates  of the stimulating electrodes, example for tripolar case=[[xe0,ye0,ze0], [xe1,ye1,ze1], [xe2,ye2,ze2]] (order is important with the middle being the cathode), INTRA case only use the position along x for IClamp
                    'amplitude': amplitudes[i], # Pulse amplitude (nA)
                    'freq': frequency[i], # Frequency of the sin pulse (kHz)
                    'duty_cycle': duty[i], # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
                    'stim_dur' : 10, # Stimulus duration (ms)
                    'waveform': waveform[i], # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
                }
                recordingParameters = {
                    "number_contact_points": 8, #Number of points on the circle constituing the cuff electrode
                    'recording_elec_pos': [recording_pos[w]], #Position of the recording electrode along axon in um, in "BIPOLAR" case the position along axons should be given as a couple [x1,x2]
                    'number_elecs': 1, #number of electrodes along the bundle
                    'dur': h.tstop, # Simulation duration (ms)
                    'rec_CAP': True, #If false means we avoid spending time using LFPy functions
                }

                unmyelinatedParameters = {
                    'name': "unmyelinated_axon", # axon name (for neuron)
                    'L': 10000, #Axon length (micrometer)
                    'diam': 1.5, #Axon diameter (micrometer)
                    'cm' : 1.0, #Specific membrane capacitance (microfarad/cm2)
                    'Ra': 200.0, #Specific axial resistance (Ohm cm)
                    'layout3D': "DEFINE_SHAPE", # either "DEFINE_SHAPE" or "PT3D" using hoc corresponding function
                    'rec_v': True, # set voltage recorders True or False
                }
                myelinatedParametersA = {
                    'name': "myelinated_axonA", # axon name (for neuron)
                    'Nnodes': Nnodes, #Number of nodes
                    'fiberD': fiberD_A, #myelinatedDistribution, Diameter of the fiber
                    'layout3D': "DEFINE_SHAPE", # either "DEFINE_SHAPE" or "PT3D" using hoc corresponding function
                    'rec_v': True, # set voltage recorders True or False
                    'nodelength' : nodelength, #set node length (um)
                    'paralength1': paralength1, #set the length of the nearest paranode segment to the node
                    'paralength2': paralength2,  #set the length of the second paranode segment followed by the internodes segments
                    'interlength': interlength, #set the length of the internode part comprising the 6 segments between two paranodes2
                }


                bundleParameters = {         #parameters for Bundle classe
                        'radius_bundle': 150.0, #um Radius of the bundle (typically 0.5-1.5mm)
                        'draw_distribution': True, #Boolean stating the distribution of fibre should be drawn
                        'number_of_axons': 1, # Number of axons in the bundle
                        'p_A': p_A[k], # Percentage of myelinated fiber type A
                        'p_C': 1-p_A[k], #Percentage of unmyelinated fiber type C
                        'myelinated_A': myelinatedParametersA, #parameters for fiber type A
                        'unmyelinated': unmyelinatedParameters, #parameters for fiber type C
                }

                Parameters1 = dict(bundleParameters, **stimulusParameters)
                Parameters = dict(Parameters1, **recordingParameters)
                
                bundle = Bundle(**Parameters)

                if recordingParameters['rec_CAP'] == False:
                    # When saving voltage to file limit the number of axons to 10. If 100 unmyelinated it produces a 1Go file, and if 100 myelinated 2Go.
                    directory = "FOR_PAPER/Voltage/myelinated/"
                    if not os.path.exists(directory):
                        os.makedirs(directory)
                    save_voltage_tofile(bundle,Parameters,directory)

                else:
                    directory = "FOR_PAPER/CAP/myelinated/"
                    if not os.path.exists(directory):
                        os.makedirs(directory)
                    save_CAP_tofile(bundle,Parameters,directory)

                bundle = None






    




