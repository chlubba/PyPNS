from neuron import h
h('load_file("noload.hoc")') # setup the environment without loading the GUI
import numpy as np
import matplotlib.pyplot as plt

# set some basic parameters
h.v_init = -65
h.celsius = 33 # set temperature in celsius
h.tstop = 3e1 # set simulation duration (ms)
h.dt = 0.0025 #0.0005 # set time step (ms)
h.finitialize() # initialize voltage state

# define one section, the axon
axon = h.Section(name='unmyelinated_axon')
axon.nseg = 200
axon.L = 10000
axon.diam = 1
axon.insert('hh')

# define a stimulus and connect it to the axon
stim = h.IClamp(0, axon)
stim.delay = 0
stim.dur = 1
stim.amp = 1

# record the membrane voltage
axonSegVoltages = h.List()
for seg in axon:
    vrec = h.Vector(int(h.tstop / h.dt + 1))
    vrec.record(seg._ref_v)
    axonSegVoltages.append(vrec)

# and the time
timeVector = h.Vector()
timeVector.record(h._ref_t)

# run the simulation
h.run()

# get the recorded quantities as numpy arrays
timeArray = np.array(timeVector)
vAxon = np.array(axonSegVoltages)

# plot it
plt.plot(timeArray, vAxon.T)
plt.xlabel('time [ms]')
plt.ylabel('segment membrane voltage [mV]')
plt.show()