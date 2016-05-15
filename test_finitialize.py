import neuron
from neuron import h
import neuron.gui
import numpy as np
from matplotlib import pyplot

h.v_init = -90

h.celsius = 33  # set temperature in celsius
h.tstop = 3e1  # set simulation duration (ms)
h.dt = 0.0025  # 0.0005 # set time step (ms)
h.finitialize()  # initialize voltage state

# create first axon and set properties
h('create axon0')
axon0 = h.axon0
# axon0.nseg = 20
# axon0.L = 1000
# axon0.diam = 1
# axon0.insert('hh')

# # create current clamp and connect to FIRST axon
# stim1 = h.IClamp(0, axon0)
# stim1.delay = 10
# stim1.dur = 1
# stim1.amp = 1

# set voltage recorder for first axon
vAxonVec0 = h.Vector()
vAxonVec0.record(axon0(0)._ref_v)

# record time, too
timeVector = h.Vector()
timeVector.record(h._ref_t)

# run simulation
h.run()

# convert recorded data to np.array for plotting
timeArray = np.array(timeVector)
vAxon0 = np.array(vAxonVec0)

# plot all signals
pyplot.plot(timeVector, vAxon0, label='First axon')
pyplot.xlabel('time [ms]')
pyplot.ylabel('voltage [mV]')
pyplot.title('When adding IClamps to every axon:')
pyplot.legend()
pyplot.show()
