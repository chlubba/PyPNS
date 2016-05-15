# --------------------------------------------------------------------
# 2/02
# Cameron C. McIntyre
# SIMULATION OF PNS MYELINATED AXON
#
# This model is described in detail in:
#
# McIntyre CC, Richardson AG, and Grill WM. Modeling the excitability of
# mammalian nerve fibers: influence of afterpotentials on the recovery
# cycle. Journal of Neurophysiology 87:995-1006, 2002.
#
# This model can not be used with NEURON v5.1 as errors in the
# extracellular mechanism of v5.1 exist related to xc. The original
# stimulations were run on v4.3.1. NEURON v5.2 has corrected the
# limitations in v5.1 and can be used to run this model.
# ----------------------------------------------------------------------

# load_proc("nrnmainmenu")

from math import pi
import numpy as np
import matplotlib.pyplot as plt

from neuron import h
h('load_file("noload.hoc")')


# set the globals
h.celsius=37
h.v_init=-80
h.dt=0.005
h.tstop=30

# Intracellular stimuluation parameters
istim=2
delay=1
pw=0.1

# topological parameters
axonnodes=21
paranodes1=40
paranodes2=40
axoninter=120
axontotal=221

# morphological parameters
fiberD=10.0 # choose from 5.7, 7.3, 8.7, 10.0, 11.5, 12.8, 14.0, 15.0, 16.0
paralength1=3
nodelength=1.0
space_p1=0.002
space_p2=0.004
space_i=0.004

# electrical parameters
rhoa=0.7e6 # Ohm-um//
mycm=0.1 # uF/cm2/lamella membrane//
mygm=0.001 # S/cm2/lamella membrane//



if (fiberD==5.7):
    g=0.605
    axonD=3.4
    nodeD=1.9
    paraD1=1.9
    paraD2=3.4
    deltax=500
    paralength2=35
    nl=80
# if (fiberD==7.3) {g=0.630 axonD=4.6 nodeD=2.4 paraD1=2.4 paraD2=4.6 deltax=750 paralength2=38 nl=100}
# if (fiberD==8.7) {g=0.661 axonD=5.8 nodeD=2.8 paraD1=2.8 paraD2=5.8 deltax=1000 paralength2=40 nl=110}
if (fiberD==10.0):
    g=0.690
    axonD=6.9
    nodeD=3.3
    paraD1=3.3
    paraD2=6.9
    deltax=1150
    paralength2=46
    nl=120
# if (fiberD==11.5) {g=0.700 axonD=8.1 nodeD=3.7 paraD1=3.7 paraD2=8.1 deltax=1250 paralength2=50 nl=130}
# if (fiberD==12.8) {g=0.719 axonD=9.2 nodeD=4.2 paraD1=4.2 paraD2=9.2 deltax=1350 paralength2=54 nl=135}
# if (fiberD==14.0) {g=0.739 axonD=10.4 nodeD=4.7 paraD1=4.7 paraD2=10.4 deltax=1400 paralength2=56 nl=140}
# if (fiberD==15.0) {g=0.767 axonD=11.5 nodeD=5.0 paraD1=5.0 paraD2=11.5 deltax=1450 paralength2=58 nl=145}
# if (fiberD==16.0) {g=0.791 axonD=12.7 nodeD=5.5 paraD1=5.5 paraD2=12.7 deltax=1500 paralength2=60 nl=150}
Rpn0=(rhoa*.01)/(pi*((((nodeD/2)+space_p1)**2)-((nodeD/2)**2)))
Rpn1=(rhoa*.01)/(pi*((((paraD1/2)+space_p1)**2)-((paraD1/2)**2)))
Rpn2=(rhoa*.01)/(pi*((((paraD2/2)+space_p2)**2)-((paraD2/2)**2)))
Rpx=(rhoa*.01)/(pi*((((axonD/2)+space_i)**2)-((axonD/2)**2)))
interlength=(deltax-nodelength-(2*paralength1)-(2*paralength2))/6


# objectvar stim
#
# create node[axonnodes], MYSA[paranodes1], FLUT[paranodes2], STIN[axoninter]
# access node[0]	//APD

# allseclist=h.SectionList()

nodes = []
for i in range(axonnodes):
    node = h.Section()
    node.nseg=1
    node.diam=nodeD
    node.L=nodelength
    node.Ra=rhoa/10000
    node.cm=2
    node.insert('axnode')

    node.insert('extracellular')
    node.xraxial[0] = Rpn0
    node.xraxial[1] = Rpn0
    node.xg[0] =1e10
    node.xg[1] =1e10
    node.xc[0] =0
    node.xc[1] =0

    nodes.append(node)
    # allseclist.append(sec=node)

MYSAs = []
for i in range(paranodes1):
    MYSA = h.Section()
    MYSA.nseg=1
    MYSA.diam=fiberD
    MYSA.L=paralength1
    MYSA.Ra=rhoa*(1/(paraD1/fiberD)**2)/10000
    MYSA.cm=2*paraD1/fiberD
    MYSA.insert('pas')
    MYSA.g_pas=0.001*paraD1/fiberD
    MYSA.e_pas=-80

    MYSA.insert('extracellular')
    MYSA.xraxial[0] = Rpn1
    MYSA.xraxial[1] = Rpn1
    MYSA.xg[0] = mygm/(nl*2)
    MYSA.xg[1] = mygm/(nl*2)
    MYSA.xc[0] = mycm/(nl*2)
    MYSA.xc[1] = mycm/(nl*2)

    MYSAs.append(MYSA)
    # allseclist.append(sec=MYSA)

FLUTs = []
for i in range(paranodes2):
    FLUT = h.Section()
    FLUT.nseg=1
    FLUT.diam=fiberD
    FLUT.L=paralength2
    FLUT.Ra=rhoa*(1/(paraD2/fiberD)**2)/10000
    FLUT.cm=2*paraD2/fiberD
    FLUT.insert('pas')
    FLUT.g_pas=0.0001*paraD2/fiberD
    FLUT.e_pas=-80

    FLUT.insert('extracellular')
    FLUT.xraxial[0] = Rpn2
    FLUT.xraxial[1] = Rpn2
    FLUT.xg[0] = mygm/(nl*2)
    FLUT.xg[1] = mygm/(nl*2)
    FLUT.xc[0] = mycm/(nl*2)
    FLUT.xc[1] = mycm/(nl*2)

    FLUTs.append(FLUT)
    # allseclist.append(sec=FLUT)

STINs = []
for i in range(axoninter):
    STIN = h.Section()
    STIN.nseg=1
    STIN.diam=fiberD
    STIN.L=interlength
    STIN.Ra=rhoa*(1/(axonD/fiberD)**2)/10000
    STIN.cm=2*axonD/fiberD
    STIN.insert('pas')
    STIN.g_pas=0.0001*axonD/fiberD
    STIN.e_pas=-80

    STIN.insert('extracellular')
    STIN.xraxial[0] = Rpx
    STIN.xraxial[1] = Rpx
    STIN.xg[0] = mygm/(nl*2)
    STIN.xg[1] = mygm/(nl*2)
    STIN.xc[0] = mycm/(nl*2)
    STIN.xc[1] = mycm/(nl*2)

    STINs.append(STIN)
    # allseclist.append(sec=STIN)

for i in range(axonnodes-1):
    MYSAs[2*i].connect(nodes[i],1,0)
    FLUTs[2*i].connect(MYSAs[2*i],1,0)
    STINs[6*i].connect(FLUTs[2*i],1,0)
    STINs[6*i+1].connect(STINs[6*i],1,0)
    STINs[6*i+2].connect(STINs[6*i+1],1,0)
    STINs[6*i+3].connect(STINs[6*i+2],1,0)
    STINs[6*i+4].connect(STINs[6*i+3],1,0)
    STINs[6*i+5].connect(STINs[6*i+4],1,0)
    FLUTs[2*i+1].connect(STINs[6*i+5],1,0)
    MYSAs[2*i+1].connect(FLUTs[2*i+1],1,0)
    nodes[i+1].connect(MYSAs[2*i+1],1,0)

# h.define_shape()

h.finitialize()
h.fcurrent()

# create current clamp and connect to FIRST axon
stim1 = h.IClamp(0, sec=nodes[10]) # sec=
stim1.delay = 15
stim1.dur = 1
stim1.amp = 1

# initialize()
#
# intracellular stimulus//

# h('node[10]{stim=new IClamp()}')
# h('allseclist[0]{stim=new IClamp()}')
# 		stim=new IClamp()
# 		stim.loc(.5)
# 		stim.del=delay
# 		stim.dur=pw
# 		stim.amp=istim
# 		}
# }
# stimul()
#
# xpanel("Stimulus parameters")
# 	xvalue("Stimulus Amplitude (nA)", "istim", 1, "stimul()", 1)
# 	xvalue("Pulse Duration (ms)", "pw", 1)
# 	xvalue("Onset Delay (ms)", "delay", 1)
# xpanel(100,100)

# set voltage recorder for first axon
vVec0 = h.Vector()
vVec0.record(nodes[10](0)._ref_v)

# record time, too
tVec = h.Vector()
tVec.record(h._ref_t)

# run simulation
h.run()

# convert recorded data to np.array for plotting
timeArray = np.array(tVec)
vAxon0 = np.array(vVec0)

# plot all signals
plt.plot(tVec, vAxon0, label='Voltage')
plt.xlabel('time [ms]')
plt.ylabel('voltage [mV]')
plt.title('When adding IClamps to every axon:')
plt.legend()
plt.show()