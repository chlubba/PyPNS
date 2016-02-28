from neuron import h
import neuron.gui
import numpy as np
from matplotlib import pyplot

h.tstop = 3e1 # set simulation duration (ms)
h.dt = 0.0025 #0.0005 # set time step (ms)
h.finitialize(-65) # initialize voltage state

axon = h.Section(name='unmyelinated_axon')
axon.nseg = 20
axon.L = 1000
axon.diam = 1

# axon.insert('extracellular')
# axon.insert('xtra')

axon.insert('hh') # insert a Hodgkin & Huxley channel

# axon.gnabar_hh = 0.120*2
# axon.gkbar_hh = 0.036
# axon.gl_hh = 0.0003
# axon.ena = 50
# axon.ek = -77
# axon.el_hh = -54.3

expSyn = h.ExpSyn(0, axon)
expSyn.e = 10
expSyn.i = 0.2
expSyn.tau = 3

vecStim = h.VecStim()
vec = h.Vector([5, 10, 15, 20, 25])#[1, 2])
vecStim.play(vec)

netCon = h.NetCon(vecStim, expSyn)
netCon.weight[0] = 1 #

vec = {}
for var in 'v_axon0','v_axon1','t', 'i_syn': #
    vec[var] = h.Vector()

vec['v_axon0'].record(axon(0)._ref_v)
vec['v_axon1'].record(axon(1)._ref_v)
vec['i_syn'].record(expSyn._ref_i)
vec['t'].record(h._ref_t)

print np.shape(range(0, int(h.tstop/h.dt)))

# for t in range(0, int(h.tstop/h.dt)):
#     #do some stuff...
#     h.fadvance()
h.run()

vAxon0 = np.array(vec['v_axon0'])
vAxon1 = np.array(vec['v_axon1'])
iSyn = np.array(vec['i_syn'])
timeArray = np.array(vec['t'])

# print np.shape(vSoma)
# print np.shape(timeArray)


pyplot.figure(figsize=(8,4)) # Default figsize is (8,6)
# pyplot.plot(timeArray, vAxon0)
# pyplot.plot(timeArray, vAxon1)
pyplot.plot(timeArray, iSyn)
pyplot.xlabel('time (ms)')
pyplot.ylabel('mV')
pyplot.show()

