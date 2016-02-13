from neuron import h
import neuron.gui
import numpy as np
from matplotlib import pyplot

h.tstop = 3e1 # set simulation duration (ms)
h.dt = 0.0025 #0.0005 # set time step (ms)
h.finitialize(-65) # initialize voltage state

soma = h.Section(name='soma')
soma.insert('hh')

print "type(soma) =", type(soma)
print "type(soma(0.5)) =", type(soma(0.5))

expSyn = h.ExpSyn(0.5, soma)
expSyn.e = 10
expSyn.i = 0.2
expSyn.tau = 3

vecStim = h.VecStim()
vec = h.Vector([10,25])#[1, 2])
vecStim.play(vec)

netCon = h.NetCon(vecStim, expSyn)
netCon.weight[0] = 1

vec = {}
for var in 'v_soma', 't':
    vec[var] = h.Vector()

vec['v_soma'].record(soma(0.5)._ref_v)
vec['t'].record(h._ref_t)

print np.shape(range(0, int(h.tstop/h.dt)))

# for t in range(0, int(h.tstop/h.dt)):
#     #do some stuff...
#     h.fadvance()
h.run()

vSoma = np.array(vec['v_soma'])
timeArray = np.array(vec['t'])

# print np.shape(vSoma)
# print np.shape(timeArray)


pyplot.figure(figsize=(8,4)) # Default figsize is (8,6)
pyplot.plot(timeArray, vSoma)
pyplot.xlabel('time (ms)')
pyplot.ylabel('mV')
pyplot.show()

