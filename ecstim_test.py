import neuron
from neuron import h
import neuron.gui
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt


def set_nsegs_lambda_f(section, frequency=100, d_lambda=0.1):
    # Set the number of segments for section according to the
    # d_lambda-rule for a given input frequency
    #     frequency: float, frequency at whihc AC length constant is computed
    #     d_lambda: float,


    section.nseg = int((section.L / (d_lambda * h.lambda_f(frequency, sec=section)) + .9)/ 2) * 2 + 1
    print "Number of segments for unmyelinated axon via d_lambda: " + str(section.nseg)

def hh_init(section, gna=0.120, gk=0.036, gl=0.0003, ena=50, ek=-77, el=-54.3):

        # g unit S/cm2 and e unit mV
        # gna: maximum specific sodium channel conductance
        # gk: maximum specific potassium channel conductance
        # gl: maximum specific leakage conductance
        # ena: reversal potential for the sodium channel
        # ek: reversal potential for the potassium channel
        # el: reversal potential for the leakage channel

        section.insert('hh') # insert a Hodgkin & Huxley channel

        section.gnabar_hh = gna
        section.gkbar_hh = gk
        section.gl_hh = gl
        section.ena = ena
        section.ek = ek
        section.el_hh = el


def rectangularStimulusSignal(stimDur=0.1, amplitude=1, frequency=10, dutyCycle=0.5, waveform='MONOPHASIC', timeRes=.0025, delay=0, invert=False):

    t = np.arange(0, stimDur, timeRes)

    if waveform == 'MONOPHASIC':
        stimulusSignal = amplitude * 0.5 * signal.square(2 * np.pi * frequency * t, duty=dutyCycle) + amplitude * 0.5
    elif waveform == 'BIPHASIC':
        stimulusSignal = amplitude * signal.square(2 * np.pi * frequency * t, duty=dutyCycle)
    else:
        print "You didn't choose the right waveform either MONOPHASIC or BIPHASIC, it has been set to default MONOPHASIC"
        stimulusSignal = amplitude * 0.5 * signal.square(2 * np.pi * frequency * t, duty=dutyCycle) + amplitude * 0.5

    if invert:
        stimulusSignal = -stimulusSignal

    # apply delay
    stimulusSignal = np.concatenate((np.zeros(delay / timeRes), stimulusSignal))
    t = np.arange(0, stimDur + delay, timeRes)

    return t, stimulusSignal


# record the membrane potentials
def set_voltage_recorders(sectionList):
    vreclist = h.List()

    for sec in sectionList:
        for seg in sec:
            vrec = h.Vector(int(h.tstop / h.dt + 1))
            vrec.record(seg._ref_v)
            vreclist.append(vrec)

    return vreclist

# set the globals
timeRes = 0.0025
# h.celsius=37
# h.v_init=-65
h.dt=timeRes
h.tstop=30

# define sections
axonSection0 = h.Section()
# axonSection0.L = 10000
# axonSection0.diam = 1
# set_nsegs_lambda_f(axonSection0)

axonSection0.nseg = 1
axonSection0.diam = 18.8
axonSection0.L = 18.8
axonSection0.Ra = 123.0

# hh_init(axonSection0)
axonSection0.insert('hh')
# axonSection0.insert('extracellular')

# list sections
sectionList = h.SectionList()
sectionList.append(axonSection0)
#
# # define stimulation signals
# t, stimulusSignal = rectangularStimulusSignal(timeRes=timeRes, delay=5)
# tVec = h.Vector(t)
#
# numberOfSegments = 0
# for sec in sectionList:
#     numberOfSegments += sec.nseg
# factors =  np.linspace(-400, 0, numberOfSegments) # np.ones(numberOfSegments)*(100) # np.concatenate((np.zeros(numberOfSegments-1), [-100]))# [0,0,-10] # np.linspace(2.,10.,numberOfSegments)

# # connect to cells
# signalVectors = []
# segCounter = 0
# for sec in sectionList:
#     for seg in sec:
#
#         stimulusSignalSection = stimulusSignal * factors[segCounter]
#         # plt.plot(stimulusSignalSection)
#         # plt.show()
#
#         signalVectors.append(h.Vector(stimulusSignalSection))
#         signalVectors[-1].play(seg._ref_e_extracellular, tVec)
#
#         segCounter += 1

# create current clamp and connect to SECOND axon
stim0 = h.IClamp(0, sec=axonSection0)
stim0.delay = 0
stim0.dur = 30
stim0.amp = 10

# setup recording
vreclist = set_voltage_recorders(sectionList)

# start
h.finitialize()
h.run()

# get voltages from segments
voltages = np.transpose(np.array(vreclist))


for i in range(voltages.shape[1]):
    plt.plot(np.arange(0,voltages.shape[0])*timeRes, voltages[:,i], label=str(i))
plt.legend()
plt.show()




