import neuron
from neuron import h
import neuron.gui
import numpy as np
from matplotlib import pyplot

h.celsius = 33 # set temperature in celsius
h.tstop = 3e1 # set simulation duration (ms)
h.dt = 0.0025 #0.0005 # set time step (ms)
h.finitialize(-65) # initialize voltage state

axons = []

# for i in range(2):
#
#     axonTest = h.Section(name='unmyelinated_axon')
#     axonTest.nseg = 20
#     axonTest.L = 1000
#     axonTest.diam = 1
#
#     axonTest.insert('hh')
#
#     axons.append(axonTest)

# for i in range (2):
#
#     axon = axons[i]
#
#     stim = h.IClamp(0, axon)
#     stim.delay = 0
#     stim.dur = 1
#     stim.amp = 100

# for i in range(2):
#     print list(axons[0].allseg())

axonTest1 = h.Section(name='unmyelinated_axon')
axonTest1.nseg = 20
axonTest1.L = 1000
axonTest1.diam = 1
#
axonTest1.insert('hh')

# stim = h.IClamp(0, axonTest1)
# stim.delay = 0
# stim.dur = 1
# stim.amp = 100

axonTest2 = h.Section(name='unmyelinated_axon')
axonTest2.nseg = 20
axonTest2.L = 1000
axonTest2.diam = 1

axonTest2.insert('hh')

stim2 = h.IClamp(1, axonTest2)
stim2.delay = 0
stim2.dur = 1
stim2.amp = 100

print list(axonTest1.allseg())
print list(axonTest2.allseg())

vAxonVec1 = h.Vector()
vAxonVec1.record(axonTest1(0)._ref_v)

vAxonVec2 = h.Vector()
vAxonVec2.record(axonTest2(0)._ref_v)

timeVector = h.Vector()
timeVector.record(h._ref_t)

h.run()

timeArray = np.array(timeVector)
vAxon1 = np.array(vAxonVec1)
vAxon2 = np.array(vAxonVec2)

pyplot.plot(timeVector,vAxon1, label='First axon')
pyplot.plot(timeVector,vAxon2, label='Second axon')

# pyplot.plot(timeVector,iSyn)
pyplot.legend()
pyplot.show()