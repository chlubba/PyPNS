import neuron
from neuron import h
import neuron.gui
import numpy as np
from matplotlib import pyplot

h.celsius = 33 # set temperature in celsius
h.tstop = 3e1 # set simulation duration (ms)
h.dt = 0.0025 #0.0005 # set time step (ms)
h.finitialize(-65) # initialize voltage state

axonTest = h.Section(name='unmyelinated_axon')
axonTest.nseg = 20
axonTest.L = 1000
axonTest.diam = 1

axonTest.insert('hh')

axonTest.insert('extracellular')
axonTest.insert('xtra')


# # self.stim = h.ExpSyn(0, self.axon.allseclist)
# stim = h.ExpSyn(0, axonTest)
# stim.e = 10
# stim.i = 0.2
# stim.tau = 3
#
# vecStim = h.VecStim()
# vec = h.Vector([5, 25])#[1, 2])
# vecStim.play(vec)
#
# netCon = h.NetCon(vecStim, stim)
# netCon.weight[0] = 1

stim = h.IClamp(0, axonTest)
stim.delay = 0
stim.dur = 1
stim.amp = 100

# firingTimes = np.array([5, 25])
#
# time = np.linspace(0, h.tstop, (h.tstop+2*h.dt)/h.dt, endpoint=True)
# signal = np.zeros(np.shape(time))
# timesArray = np.array(firingTimes/h.dt)
# print timesArray
# signal[timesArray.astype(int)] = np.ones(np.shape(firingTimes))
#
# # pyplot.plot(time,signal)
# # pyplot.show()
#
# vec = h.Vector(signal)
#
# vec.play(stim._ref_amp,h.dt)
# # vec.play(stim.amp, h.dt)

#observe synapse current for testing purposes
i_syn_vec = h.Vector()
i_syn_vec.record(stim._ref_i)

vAxonVec = h.Vector()
vAxonVec.record(axonTest(0)._ref_v)

timeVector = h.Vector()
timeVector.record(h._ref_t)

h.run()

# just shortly test if synapse currents happen
iSyn = np.array(i_syn_vec)
timeArray = np.array(timeVector)
vAxon = np.array(vAxonVec)

pyplot.plot(timeVector,vAxon)
# pyplot.plot(timeVector,iSyn)
pyplot.show()