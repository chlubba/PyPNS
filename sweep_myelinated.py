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

# BLOCK inserted
# 'name': "myelinated_axonA", # axon name (for neuron)
# 'Nnodes': 11, #Number of nodes
# 'fiberD': fiberD_A, #myelinatedDistribution, Diameter of the fiber
# 'layout3D': "DEFINE_SHAPE", # either "DEFINE_SHAPE" or "PT3D" using hoc corresponding function
# 'rec_v': True, # set voltage recorders True or False
# 'nodelength' : nodelength, #set node length (um)
# 'paralength1': paralength1, #set the length of the nearest paranode segment to the node
# 'paralength2': paralength2_A,  #set the length of the second paranode segment followed by the internodes segments
# 'interlength': interlength_A, #set the length of the internode part comprising the 6 segments between two paranodes2

layout3D = "DEFINE_SHAPE"
axonnodes = 11
coord = [0,0]
fiberD = 5.7 #choose from 5.7, 7.3, 8.7, 10.0, 11.5, 12.8, 14.0, 15.0, 16.0

print 'Myelinated fiber diameter: ' + str(fiberD)
# topological parameters
"""axonnodes=21"""
paranodes1= 2*(axonnodes-1)
paranodes2= 2*(axonnodes-1)
axoninter= 6*(axonnodes-1)
axontotal= axonnodes+paranodes1+paranodes2+axoninter

# morphological parameters
paralength1=3
nodelength=1.0
space_p1=0.002
space_p2=0.004
space_i=0.004

#electrical parameters
rhoa=0.7e6 #Ohm-um
mycm=0.1 #uF/cm2/lamella membrane
mygm=0.001 #S/cm2/lamella membrane


###### from proc dependent_var() #############
if (fiberD==5.7):
    g=0.605
    axonD=3.4
    nodeD=1.9
    paraD1=1.9
    paraD2=3.4
    deltax=500
    paralength2=35
    nl=80
if (fiberD==7.3):
    g=0.630
    axonD=4.6
    nodeD=2.4
    paraD1=2.4
    paraD2=4.6
    deltax=750
    paralength2=38
    nl=100
if (fiberD==8.7):
    g=0.661
    axonD=5.8
    nodeD=2.8
    paraD1=2.8
    paraD2=5.8
    deltax=1000
    paralength2=40
    nl=110
if (fiberD==10.0):
    g=0.690
    axonD=6.9
    nodeD=3.3
    paraD1=3.3
    paraD2=6.9
    deltax=1150
    paralength2=46
    nl=120
if (fiberD==11.5):
    g=0.700
    axonD=8.1
    nodeD=3.7
    paraD1=3.7
    paraD2=8.1
    deltax=1250
    paralength2=50
    nl=130
if (fiberD==12.8):
    g=0.719
    axonD=9.2
    nodeD=4.2
    paraD1=4.2
    paraD2=9.2
    deltax=1350
    paralength2=54
    nl=135
if (fiberD==14.0):
    g=0.739
    axonD=10.4
    nodeD=4.7
    paraD1=4.7
    paraD2=10.4
    deltax=1400
    paralength2=56
    nl=140
if (fiberD==15.0):
    g=0.767
    axonD=11.5
    nodeD=5.0
    paraD1=5.0
    paraD2=11.5
    deltax=1450
    paralength2=58
    nl=145
if (fiberD==16.0):
    g=0.791
    axonD=12.7
    nodeD=5.5
    paraD1=5.5
    paraD2=12.7
    deltax=1500
    paralength2=60
    nl=150


Rpn0=(rhoa*0.01)/(math.pi*((math.pow((nodeD/2)+space_p1,2))-(math.pow(nodeD/2,2))))
Rpn1=(rhoa*0.01)/(math.pi*((math.pow((paraD1/2)+space_p1,2))-(math.pow(paraD1/2,2))))
Rpn2=(rhoa*0.01)/(math.pi*((math.pow((paraD2/2)+space_p2,2))-(math.pow(paraD2/2,2))))
Rpx=(rhoa*0.01)/(math.pi*((math.pow((axonD/2)+space_i,2))-(math.pow(axonD/2,2))))
interlength=(deltax-nodelength-(2*paralength1)-(2*paralength2))/6

#### from initialize() ####
nodes = []
for i in range(axonnodes):
    node = h.Section()
    node.nseg = 1
    node.diam = nodeD
    node.L = nodelength
    node.Ra = rhoa/10000
    node.cm = 2
    node.insert('axnode')
    #h('insert extracellular xraxial=Rpn0 xg=1e10 xc=0')
    node.insert('extracellular')
    node.xraxial[0] = Rpn0
    node.xraxial[1] = Rpn0
    node.xg[0] =1e10
    node.xg[1] =1e10
    node.xc[0] =0
    node.xc[1] =0
    # node.insert('xtra')

    nodes.append(node)

MYSAs = [] # paranodes1
for i in range(paranodes1):
    MYSA = h.Section()
    MYSA.nseg = 1
    MYSA.diam = fiberD
    MYSA.L = paralength1
    MYSA.Ra = rhoa*(1/math.pow(paraD1/fiberD,2))/10000
    MYSA.cm = 2*paraD1/fiberD
    MYSA.insert('pas')
    MYSA.g_pas = 0.001*paraD1/fiberD
    MYSA.e_pas = -80
    #h('insert extracellular xraxial=Rpn1 xg=mygm/(nl*2) xc=mycm/(nl*2)')
    MYSA.insert('extracellular')
    MYSA.xraxial[0] = Rpn1
    MYSA.xraxial[1] = Rpn1
    MYSA.xg[0] = mygm/(nl*2)
    MYSA.xg[1] = mygm/(nl*2)
    MYSA.xc[0] = mycm/(nl*2)
    MYSA.xc[1] = mycm/(nl*2)
    # MYSA.insert('xtra')

    MYSAs.append(MYSA)

FLUTs = [] # paranodes2
for i in range(paranodes2):
    FLUT = h.Section()
    FLUT.nseg = 1
    FLUT.diam = fiberD
    FLUT.L = paralength2
    FLUT.Ra = rhoa*(1/math.pow(paraD2/fiberD,2))/10000
    FLUT.cm = 2*paraD2/fiberD
    FLUT.insert('pas')
    FLUT.g_pas = 0.0001*paraD2/fiberD
    FLUT.e_pas = -80
    #h('insert extracellular xraxial=Rpn2 xg=mygm/(nl*2) xc=mycm/(nl*2)')
    FLUT.insert('extracellular')
    FLUT.xraxial[0] =Rpn2
    FLUT.xraxial[1] =Rpn2
    FLUT.xg[0] = mygm/(nl*2)
    FLUT.xg[1] = mygm/(nl*2)
    FLUT.xc[0] = mycm/(nl*2)
    FLUT.xc[1] = mycm/(nl*2)
    # FLUT.insert('xtra')

    FLUTs.append(FLUT)

STINs = [] # internodes
for i in range(axoninter):
    STIN = h.Section()
    STIN.nseg = 1
    STIN.diam = fiberD
    STIN.L = interlength
    STIN.Ra = rhoa*(1/math.pow(axonD/fiberD,2))/10000
    STIN.cm = 2*axonD/fiberD
    STIN.insert('pas')
    STIN.g_pas = 0.0001*axonD/fiberD
    STIN.e_pas = -80
    #h('insert extracellular xraxial=Rpx xg=mygm/(nl*2) xc=mycm/(nl*2)')
    STIN.insert('extracellular')
    STIN.xraxial[0] = Rpx
    STIN.xraxial[1] = Rpx
    STIN.xg[0] = mygm/(nl*2)
    STIN.xg[1] = mygm/(nl*2)
    STIN.xc[0] = mycm/(nl*2)
    STIN.xc[1] = mycm/(nl*2)
    # STIN.insert('xtra')

    STINs.append(STIN)
if (layout3D == "DEFINE_SHAPE"):
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
    h.define_shape()
elif (layout3D == "PT3D"):
    ### PLACE AXON SPATIALLY #############################################################################################################################################################
    for i in range(axonnodes-1):
        h.pt3dclear(sec=nodes[i])
        h.pt3dadd(0, coord[0], coord[1], nodes[i].diam, sec=nodes[i])
        h.pt3dadd(nodelength*i+interlength*6*i+2*i*paralength1+ 2*i*paralength2 + nodelength, coord[0], coord[1], nodes[i].diam, sec=nodes[i])

        MYSAs[2*i].connect(nodes[i],1,0)
        h.pt3dclear(sec=MYSAs[2*i])
        h.pt3dadd(nodelength*i+interlength*6*i+2*i*paralength1+ 2*i*paralength2 + nodelength, coord[0], coord[1], MYSAs[2*i].diam, sec=MYSAs[2*i])
        h.pt3dadd(nodelength*(i+1)+interlength*6*i+2*i*paralength1+  2*i*paralength2 +paralength1, coord[0], coord[1], MYSAs[2*i].diam, sec=MYSAs[2*i])

        FLUTs[2*i].connect(MYSAs[2*i],1,0)
        h.pt3dclear(sec=FLUTs[2*i])
        h.pt3dadd(nodelength*(i+1)+interlength*6*i+2*i*paralength1+  2*i*paralength2 +paralength1, coord[0], coord[1], FLUTs[2*i].diam, sec=FLUTs[2*i])
        h.pt3dadd(nodelength*(i+1)+interlength*6*i+2*i*paralength1+  2*i*paralength2 + paralength1 + paralength2, coord[0], coord[1], FLUTs[2*i].diam, sec=FLUTs[2*i])
        for j in range(6):
            STINs[6*i+j].connect(FLUTs[2*i],1,0)
            h.pt3dclear(sec=STINs[6*i+j])
            h.pt3dadd(nodelength*(i+1)+interlength*6*i+2*i*paralength1+  2*i*paralength2 + paralength1 + paralength2+ j*interlength, coord[0], coord[1], STINs[6*i+j].diam, sec=STINs[6*i+j])
            h.pt3dadd(nodelength*(i+1)+interlength*6*i+2*(i+1.0/2)*paralength1+  2*(i+1.0/2)*paralength2 + (j+1)*interlength, coord[0], coord[1], STINs[6*i+j].diam, sec=STINs[6*i+j])
        FLUTs[2*i+1].connect(STINs[6*i+5],1,0)
        h.pt3dclear(sec=FLUTs[2*i+1])
        h.pt3dadd(nodelength*(i+1)+interlength*6*i+2*(i+1.0/2)*paralength1+  2*(i+1.0/2)*paralength2 + 6*interlength, coord[0], coord[1], FLUTs[2*i+1].diam, sec=FLUTs[2*i+1])
        h.pt3dadd(nodelength*(i+1)+interlength*6*(i+1)+2*(i+1.0/2)*paralength1+  2*(i+1.0/2)*paralength2 +paralength1, coord[0], coord[1], FLUTs[2*i+1].diam, sec=FLUTs[2*i+1])
        MYSAs[2*i+1].connect(FLUTs[2*i+1],1,0)
        h.pt3dclear(sec=MYSAs[2*i+1])
        h.pt3dadd(nodelength*(i+1)+interlength*6*(i+1)+2*(i+1.0/2)*paralength1+  2*(i+1.0/2)*paralength2 +paralength1, coord[0], coord[1], MYSAs[2*i+1].diam, sec=MYSAs[2*i+1])
        h.pt3dadd(nodelength*(i+1)+interlength*6*(i+1)+2*(i+1.0/2)*paralength1+  2*(i+1.0/2)*paralength2 + paralength1 + paralength2, coord[0], coord[1], MYSAs[2*i+1].diam, sec=MYSAs[2*i+1])
        nodes[i+1].connect(MYSAs[2*i+1],1,0)

    h.pt3dclear(sec=nodes[axonnodes-1])
    h.pt3dadd(nodelength*(axonnodes-1)+interlength*(axonnodes-1)*6+2*(axonnodes-1)*paralength1+ 2*(axonnodes-1)*paralength2, coord[0], coord[1], nodes[axonnodes-1].diam, sec=nodes[axonnodes-1])
    h.pt3dadd(nodelength*(axonnodes-1)+interlength*(axonnodes-1)*6+2*(axonnodes-1)*paralength1+ 2*(axonnodes-1)*paralength2 +nodelength, coord[0], coord[1], nodes[axonnodes-1].diam, sec=nodes[axonnodes-1])
#######################################################################################################################################################################################################
else:
    raise NameError('layout3D only "DEFINE_SHAPE" or "PT3D"')

# super(Myelinated,self).__init__(layout3D, rec_v)
# totnsegs = calc_totnsegs()
# for sec_id in allseclist:
#     for seg in sec_id:
#         h.setpointer(seg._ref_i_membrane, 'im', seg.xtra)
#         h.setpointer(seg._ref_e_extracellular, 'ex', seg.xtra)
# interpxyz()
# collect_geometry()
# BLOCK inserted


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

