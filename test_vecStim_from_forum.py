import neuron
from neuron import h
import neuron.gui
import numpy as np
from matplotlib import pyplot
#h.load_file('noload.hoc')

# set neuron simulation parameters
h.celsius = 33 # set temperature in celsius
h.tstop = 3e1 # set simulation duration (ms)
h.dt = 0.0025 #0.0005 # set time step (ms)
#h.finitialize(-100) #65 # initialize voltage state


# define axon
axonTest = h.Section(name='unmyelinated_axon')
axonTest.nseg = 20
axonTest.L = 1000
axonTest.diam = 1


# insert channels. Hodgkin-Huxley needed
axonTest.insert('hh')

# here the problems arise

# extracellular
# axonTest.insert('extracellular')

# # extracellular stimulation
# axonTest.insert('xtra')
# tsvec = h.Vector(6)
# isvec = h.Vector(6)
# tsvec.x[0] = 0.0
# tsvec.x[1] = 1.0
# tsvec.x[2] = 1.0
# tsvec.x[3] = 2.0
# tsvec.x[4] = 2.0
# tsvec.x[5] = 3.0
# isvec.x[0] = 0.0
# isvec.x[1] = 0.0
# isvec.x[2] = 1.0
# isvec.x[3] = 1.0
# isvec.x[4] = 0.0
# isvec.x[5] = 0.0
#
# # xtra's ex must be coupled to e_extracellular in the same segment
# # and its im must be coupled to i_membrane in the same segment
# for seg in axonTest:
#   # couple ex_xtra to e_extracellular
#   # hoc syntax:  setpointer ex_xtra(x), e_extracellular(x)
#   h.setpointer(axonTest(seg.x)._ref_e_extracellular, 'ex', axonTest(seg.x).xtra)
#   # couple im_xtra to i_membrane
#   # hoc syntax:  setpointer im_xtra(x), i_membrane(x)
#   h.setpointer(axonTest(seg.x)._ref_i_membrane, 'im', axonTest(seg.x).xtra)

# configure ExpSyn synapse
stim = h.ExpSyn(0, axonTest)
stim.e = 10
stim.i = 0.2
stim.tau = 3

# configure input to synapse
vecStim = h.VecStim()
vec = h.Vector([5, 25])#[1, 2])
vecStim.play(vec)

# connect synapse and VecStim input
netCon = h.NetCon(vecStim, stim)
netCon.weight[0] = 1

#observe synapse current
iSynVec = h.Vector()
iSynVec.record(stim._ref_i)

# and observe axon potentials
vAxonVec0 = h.Vector()
vAxonVec0.record(axonTest(0)._ref_v)
vAxonVec1 = h.Vector()
vAxonVec1.record(axonTest(1)._ref_v)

# record time, too
timeVector = h.Vector()
timeVector.record(h._ref_t)

# run simulation
h.finitialize() #65 # initialize voltage state
h.run()

# convert recorded signals to numpy arrays
timeArray = np.array(timeVector)
iSynArr = np.array(iSynVec)
vAxonArr0 = np.array(vAxonVec0)
vAxonArr1 = np.array(vAxonVec1)

# plot synapse currents and axon potentials
pyplot.plot(timeVector,vAxonArr0, label='Axon potential at start-position')
pyplot.plot(timeVector,vAxonArr1, label='Axon potential at end-position')
pyplot.plot(timeVector,iSynArr, label='Synapse current')
pyplot.legend()
pyplot.show()