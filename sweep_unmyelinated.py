from neuron import h
import neuron.gui
import numpy as np
from matplotlib import pyplot
from math import floor

h.tstop = 3e1 # set simulation duration (ms)
h.dt = 0.0025 #0.0005 # set time step (ms)
h.finitialize(-65) # initialize voltage state

# soma = h.Section(name='soma')
# soma.insert('hh')

axon = h.Section(name='unmyelinated_axon')
axon.nseg = 20
axon.L = 1000
axon.diam = 1

# axon.insert('extracellular')
# axon.insert('xtra')

axon.insert('hh') # insert a Hodgkin & Huxley channel

expSyn = h.ExpSyn(0, axon)
expSyn.e = 10
expSyn.i = 100#0.2
expSyn.tau = 3

vecStim = h.VecStim()
vec = h.Vector([10,25])#[1, 2])
vecStim.play(vec)

netCon = h.NetCon(vecStim, expSyn)
netCon.weight[0] = 1 #


gnabar_hh_0 = 0.120
gkbar_hh_0 = 0.036
gl_hh_0 = 0.0003
ena_0 = 50
ek_0 = -77
el_hh_0 = -54.3

# Two subplots, the axes array is 1-d
f, axarr = pyplot.subplots(6,3) #, sharex=True

for i in range(0,18,1):

    try:
        del vec
    except:
        pass

    vec = {}
    for var in 'v_axon0','v_axon1','t', 'i_syn': #
        vec[var] = h.Vector()

    vec['v_axon0'].record(axon(0)._ref_v)
    vec['v_axon1'].record(axon(1)._ref_v)
    vec['i_syn'].record(expSyn._ref_i)
    vec['t'].record(h._ref_t)

    step = float(i-5)/100#float(i-5)/100

    axon.gnabar_hh = gnabar_hh_0#*(1+step*50)
    axon.gkbar_hh = gkbar_hh_0#*(1+step*19)
    axon.gl_hh = gl_hh_0#*(1+step*19)
    axon.ena = ena_0#*(1+step*19)
    axon.ek = ek_0#*(1+step*19)
    axon.el_hh = el_hh_0*(1+step*19)

    h.run()

    vAxon0 = np.array(vec['v_axon0'])
    vAxon1 = np.array(vec['v_axon1'])
    timeArray = np.array(vec['t'])


    #axarr[i].set_title('Sharing X axis')

    axarr[floor(i/3),i%3].plot(timeArray, vAxon0, label='position 0')
    axarr[floor(i/3),i%3].plot(timeArray, vAxon1, label='position 1')
    axarr[floor(i/3),i%3].set_title('el_hh = ' + str(axon.el_hh))


    if i == 0:
        axarr[floor(i/3), i%3].legend()

    #axarr[i].xlabel('time (ms)')
    #axarr[i].ylabel('mV')

    #pyplot.figure(figsize=(8,4)) # Default figsize is (8,6)


pyplot.show()

