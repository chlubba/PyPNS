from neuron import h
# import neuron.gui
h('load_file("noload.hoc")')
import numpy as np
import math
from scipy import signal
import matplotlib.pyplot as plt

# Equivalent methods interpxyz and setrx from the xtra mechanism available on the NEURON website from Ted Carnevale
# Setrx has been modified to integrate the use of multipolar electrodes
def interpxyz(axon):
    # interpolated data, spaced at regular intervals

    # First, need to interpolate centers unto all compartments; from interpxyz.hoc
    for sec in h.allsec():

        if h.ismembrane('xtra',sec=sec):

            nn = int(h.n3d(sec=sec))
            xx = h.Vector(nn)
            yy = h.Vector(nn)
            zz = h.Vector(nn)
            length = h.Vector(nn)

            for ii in xrange(nn):
                xx.x[ii] = h.x3d(ii,sec=sec)
                yy.x[ii] = h.y3d(ii,sec=sec)
                zz.x[ii] = h.z3d(ii,sec=sec)
                length.x[ii] = h.arc3d(ii,sec=sec)

            # to use Vector class's .interpolate() must first scale the independent variable i.e. normalize length along centroid
            length.div(length.x[nn-1])
            # initialize the destination "independent" vector
            rr = h.Vector(sec.nseg+2)
            rr.indgen(1./sec.nseg)
            rr.sub(1./(2.*sec.nseg))
            rr.x[0]=0.
            rr.x[sec.nseg+1]=1.

            # length contains the normalized distances of the pt3d points along the centroid of the section.
            # These are spaced at irregular intervals.
            # range contains the normalized distances of the nodes along the centroid of the section.
            # These are spaced at regular intervals.

            # Ready to interpolate.
            xint = h.Vector(sec.nseg+2)
            yint = h.Vector(sec.nseg+2)
            zint = h.Vector(sec.nseg+2)
            xint.interpolate(rr, length, xx)
            yint.interpolate(rr, length, yy)
            zint.interpolate(rr, length, zz)

            # for each node, assign the xyz values to x_xtra, y_xtra, z_xtra
            # don't bother computing coords of the 0 and 1 ends
            # also avoid writing coords of the 1 end into the last internal node's coords
            for ii in range(1,sec.nseg+1):
                xr = rr.x[ii]
                sec(xr).x_xtra = xint.x[ii]
                sec(xr).y_xtra = yint.x[ii]
                sec(xr).z_xtra = zint.x[ii]

def setrx(axon, stim_elec, axon_pos):
    stimulation_mode = len(stim_elec) #1: monopolar,2: bipolar,3: tripolar, above just considered as multiple monopolar
    r = np.zeros(stimulation_mode)
    # now expects xyc coords as arguments

    for sec in h.allsec():

        if h.ismembrane('xtra',sec=sec):

            for seg in sec:

                for j in range(stimulation_mode):

                    [xe,ye,ze] = stim_elec[j]

                    #avoid nodes at 0 and 1 ends, so as not to override values at internal nodes
                    r[j] = math.sqrt(math.pow(seg.x_xtra - xe,2) + math.pow(seg.y_xtra-axon_pos[0] - ye,2) + math.pow(seg.z_xtra-axon_pos[1] - ze,2))

                    # 0.01 converts rho's cm to um and ohm to megohm
                    # if electrode is exactly at a node, r will be 0
                    # this would be meaningless since the location would be inside the cell
                    # so force r to be at least as big as local radius

                    if (r[j]==0):
                        r[j] = seg.diam/2.0
                if stimulation_mode == 1:
                    seg.xtra.rx = (sec.Ra / 4.0 / math.pi)*(1/r)*0.01
                elif stimulation_mode == 2:
                    seg.xtra.rx = (sec.Ra / 4.0 / math.pi)*(1/r[0]-1/r[1])*0.01
                elif stimulation_mode == 3:
                    seg.xtra.rx = (sec.Ra / 4.0 / math.pi)*(1/r[0]-1/r[1]+1/r[2])*0.01
                else:
                    sum_r = 0
                    for j in range(stimulation_mode):
                        sum_r += 1/r[j]
                    seg.xtra.rx = (sec.Ra / 4.0 / math.pi)*(1/sum_r)*0.01

# set simulation params
tStop=50.
timeRes=0.005

# stimulation params
stimDur = 1.
amplitude = 0.1 # 0.01
frequency = 1

t = np.arange(0, stimDur, timeRes)
stimSignal = amplitude * 0.5 * signal.square(2 * np.pi * frequency * t) # + amplitude * 0.5
stimSignal = np.concatenate((stimSignal, [0]))
svec = h.Vector(stimSignal)

#set up NEURON
h.v_init = -65
h.celsius = 33  # set temperature in celsius
h.tstop = tStop  # set simulation duration (ms)
h.dt = timeRes  # 0.0005 # set time step (ms)
h.finitialize()  # initialize voltage state

# create first axon and set properties
h('create axon0')
axon0 = h.axon0
axon0.insert('extracellular')
axon0.insert('xtra')
axon0.nseg = 21
axon0.L = 1000
axon0.diam = 1
axon0.insert('hh')
h.define_shape()

for seg in axon0:
    h.setpointer(seg._ref_i_membrane, 'im', seg.xtra)
    h.setpointer(seg._ref_e_extracellular, 'ex', seg.xtra)

interpxyz(axon0)
setrx(axon0, [[0, 1, 1]], [0, 0])

# set voltage recorder for first axon
vAxonVec0 = h.Vector()
vAxonVec0.record(axon0(0)._ref_v)

# record time, too
timeVector = h.Vector()
timeVector.record(h._ref_t)

# set extracellular stimulation
svec.play(h._ref_is_xtra, timeRes)

# run simulation
h.run()

# convert recorded data to np.array for plotting
timeArray = np.array(timeVector)
vAxon0 = np.array(vAxonVec0)

# plot all signals
plt.plot(timeVector, vAxon0, label='First axon')
# plt.plot(t, stimSignal, label='Stimulation signal')
plt.xlabel('time [ms]')
plt.ylabel('voltage [mV]')
plt.title('Constant voltage on stimulating electrode')
plt.legend()
plt.show()