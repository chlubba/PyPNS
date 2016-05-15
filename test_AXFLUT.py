import neuron
from neuron import h
import neuron.gui
import numpy as np
from matplotlib import pyplot

h.celsius = 33  # set temperature in celsius
h.tstop = 3e1  # set simulation duration (ms)
h.dt = 0.0025  # 0.0005 # set time step (ms)
h.finitialize(-65)  # initialize voltage state

# create first axon and set properties
h('create axon0')
axon0 = h.axon0
axon0.nseg = 20
axon0.L = 1000
axon0.diam = 1
axon0.insert('axflut') # ('axnode') #
# axon0.insert('pas')

# create current clamp and connect to FIRST axon
stim1 = h.IClamp(0, axon0)
stim1.delay = 0
stim1.dur = 1
stim1.amp = 1

# set voltage recorder for first axon
vAxonVec0 = h.Vector()
vAxonVec0.record(axon0(0)._ref_v)

# record time, too
timeVector = h.Vector()
timeVector.record(h._ref_t)

# # record currents
# iAxonVec0 = h.Vector()
# iAxonVec0.record(axon0(0)._ref_imem)

h('forall for (x,0) if (ismembrane("axflut")) gkbar_axflut(x) = 50')

# for sec in neuron.h.allsec():
#     if neuron.h.ismembrane('axflut'):
#         sec.gkbar = 0.5

# run simulation
h.run()

# convert recorded data to np.array for plotting
timeArray = np.array(timeVector)
vAxon0 = np.array(vAxonVec0)

# plot all signals
pyplot.plot(timeVector, vAxon0, label='Voltage')
# pyplot.plot(timeVector, iAxon0, label='Current')
pyplot.xlabel('time [ms]')
pyplot.ylabel('voltage [mV]')
pyplot.title('When adding IClamps to every axon:')
pyplot.legend()
pyplot.show()
