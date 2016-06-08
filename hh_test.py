from neuron import h
h('load_file("noload.hoc")')
import numpy as np
import matplotlib.pyplot as plt

h.tstop = 300

soma = h.Section(name='soma')

soma.nseg = 1
soma.diam = 18.8
soma.L = 18.8
soma.Ra = 123.0

soma.insert('hh')

# objectvar stim
# soma stim = new IClamp(0.5)
#
# stim.del = 100
# stim.dur = 100
# stim.amp = 0.1


vreclist = h.List()
for seg in soma:
    vrec = h.Vector(int(h.tstop / h.dt + 1))
    vrec.record(seg._ref_v)
    vreclist.append(vrec)


# create current clamp and connect to SECOND axon
stim0 = h.IClamp(0, sec=soma)
stim0.delay = 100
stim0.dur = 100
stim0.amp = 0.1

h.finitialize()
h.run()

# get voltages from segments
voltage = np.transpose(np.array(vreclist))

plt.plot(voltage)
plt.show()
