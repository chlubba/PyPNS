from neuron import h
h('load_file("noload.hoc")')
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

def set_nsegs_lambda_f(section, frequency=100, d_lambda=0.1):
    # Set the number of segments for section according to the
    # d_lambda-rule for a given input frequency
    #     frequency: float, frequency at whihc AC length constant is computed
    #     d_lambda: float,

    section.nseg = int((section.L / (d_lambda * h.lambda_f(frequency, sec=section)) + .9)/ 2) * 2 + 1
    print "Number of segments for unmyelinated axon via d_lambda: " + str(section.nseg)

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

timeRes = h.dt
h.tstop = 300


axon = h.Section(name='axon')
axon.diam = 1 # 18.8
axon.L = 1000 # 18.8
axon.Ra = 123.0
set_nsegs_lambda_f(axon)
axon.insert('extracellular')
# axon.insert('xtra')
axon.insert('hh')

vreclist = h.List()
for seg in axon:
    vrec = h.Vector(int(h.tstop / h.dt + 1))
    vrec.record(seg._ref_v)
    vreclist.append(vrec)


# # create current clamp and connect to axon
# stim0 = h.IClamp(0, sec=axon)
# stim0.delay = 5
# stim0.dur = 0.1
# stim0.amp = 5


# define stimulation signals
t, stimulusSignal = rectangularStimulusSignal(timeRes=timeRes, delay=5)
tVec = h.Vector(t)

sectionList = h.SectionList()
sectionList.append(axon)

numberOfSegments = 0
for sec in sectionList:
    numberOfSegments += sec.nseg
# factors =  np.linspace(-4000, 0, numberOfSegments) # np.ones(numberOfSegments)*(100) # np.concatenate((np.zeros(numberOfSegments-1), [-100]))# [0,0,-10] # np.linspace(2.,10.,numberOfSegments)
# offset = -3600
# factors =  np.concatenate((np.linspace(-4000+offset, offset, 100), np.ones(900)*offset))
# factors = np.concatenate(([-200], np.zeros(numberOfSegments-1)))
# factors = np.ones(numberOfSegments)*(-100000)
endSegment = 10
# amp = -4000/(41-endSegment)
# factors =  np.concatenate((np.linspace(-1000, 0, endSegment), np.zeros(numberOfSegments - endSegment)))
factors = np.ones(numberOfSegments)*(-10)
factors = np.cumsum(factors)*axon.L/axon.nseg

# connect to cells
signalVectors = []
segCounter = 0
for sec in sectionList:
    for seg in sec:

        stimulusSignalSection = stimulusSignal * factors[segCounter]
        # plt.plot(stimulusSignalSection)
        # plt.show()

        signalVectors.append(h.Vector(stimulusSignalSection))
        signalVectors[-1].play(seg._ref_e_extracellular, tVec)

        segCounter += 1

h.finitialize()
h.run()

# get voltages from segments
voltage = np.transpose(np.array(vreclist))

# plt.plot(voltage[:,np.linspace(0,voltage.shape[1]-1,10).astype(int)])
# plt.show()

for i in range(0,voltage.shape[1],5):
    plt.plot(np.arange(0,voltage.shape[0])*timeRes, voltage[:,i], label=str(i))
plt.legend()
plt.show()
