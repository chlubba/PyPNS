import neuron
from neuron import h
h('load_file("noload.hoc")')
import run_simulation
import refextelectrode
import numpy as np # for arrays managing
import math
import time
import glob
import os
import shutil
import copy

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

from nameSetters import getDirectoryName
"""
    Some methods in the Axon class are based on existing methods in the Python package LFPY
"""
class Axon(object):
    """
    Own initialization method of the superclass Axon
    rec_v: set voltage recorders True or False
    layout3D: either "DEFINE_SHAPE" or "PT3D" using hoc corresponding function
    # "PT3D" option has not been tested
    """

    def __init__(self,layout3D,rec_v):
        self.layout3D = layout3D
        self.rec_v = rec_v
        # LFPy initilizations        
        self.verbose = False
        self.dotprodresults = None # Not in class cell of LFPY but present in run_simulation as cell.dotprodresults
        # self.create_sectionlists()
        # self.allseclist = []
        #self.totnsegs = self.calc_totnsegs()
        
    # def create_sectionlists(self):
    #     '''Create section lists for different kinds of sections'''
    #     #list with all sections
    #     self.allsecnames = []
    #     self.allseclist = self.sections
    #     # self.allseclist = h.SectionList()
    #     # for sec in neuron.h.allsec():
    #     #     self.allsecnames.append(sec.name())
    #     #     self.allseclist.append(sec=sec)

    def calc_totnsegs(self):
        '''Calculate the number of segments in the allseclist'''
        i = 0
        for sec in self.allseclist:
            i += sec.nseg

        return i
    
    def collect_geometry(self):
        '''Collects x, y, z-coordinates from NEURON'''
        #None-type some attributes if they do not exis:
        if not hasattr(self, 'xstart'):
            self.xstart = None
            self.ystart = None
            self.zstart = None
            self.xend = None
            self.yend = None
            self.zend = None
            self.area = None
            self.diam = None
            self.length = None

        self.collect_geometry_neuron()
        self.calc_midpoints()
        
    def collect_geometry_neuron(self):
        '''Loop over allseclist to determine area, diam, xyz-start- and
        endpoints, embed geometry to self object'''
        areavec = np.zeros(self.totnsegs)
        diamvec = np.zeros(self.totnsegs)
        lengthvec = np.zeros(self.totnsegs)
        xstartvec = np.zeros(self.totnsegs)
        xendvec = np.zeros(self.totnsegs)
        ystartvec = np.zeros(self.totnsegs)
        yendvec = np.zeros(self.totnsegs)
        zstartvec = np.zeros(self.totnsegs)
        zendvec = np.zeros(self.totnsegs)

        counter = 0

        #loop over all segments
        for sec in self.allseclist:
            n3d = int(h.n3d())
            nseg = sec.nseg
            gsen2 = 1./2/nseg
            if n3d > 0:
                #create interpolation objects for the xyz pt3d info:
                L = np.zeros(n3d)
                x = np.zeros(n3d)
                y = np.zeros(n3d)
                z = np.zeros(n3d)
                for i in range(n3d):
                    L[i] = h.arc3d(i)
                    x[i] = h.x3d(i)
                    y[i] = h.y3d(i)
                    z[i] = h.z3d(i)
                #normalize as seg.x [0, 1]
                L /= sec.L

                #temporary store position of segment midpoints
                segx = np.zeros(nseg)
                for i, seg in enumerate(sec):
                    segx[i] = seg.x

                #can't be >0 which may happen due to NEURON->Python float transfer:
                segx0 = (segx - gsen2).round(decimals=6)
                segx1 = (segx + gsen2).round(decimals=6)

                #fill vectors with interpolated coordinates of start and end points
                xstartvec[counter:counter+nseg] = np.interp(segx0, L, x)
                xendvec[counter:counter+nseg] = np.interp(segx1, L, x)

                ystartvec[counter:counter+nseg] = np.interp(segx0, L, y)
                yendvec[counter:counter+nseg] = np.interp(segx1, L, y)

                zstartvec[counter:counter+nseg] = np.interp(segx0, L, z)
                zendvec[counter:counter+nseg] = np.interp(segx1, L, z)

                #fill in values area, diam, length
                for i, seg in enumerate(sec):
                    areavec[counter] = h.area(seg.x)
                    diamvec[counter] = seg.diam
                    lengthvec[counter] = sec.L/nseg
                    counter += 1
        #set self attributes
        self.xstart = xstartvec
        self.ystart = ystartvec
        self.zstart = zstartvec
        self.xend = xendvec
        self.yend = yendvec
        self.zend = zendvec
        self.area = areavec
        self.diam = diamvec
        self.length = lengthvec

    def calc_midpoints(self):
        '''Calculate midpoints of each segment'''
        self.xmid = .5*(self.xstart+self.xend).flatten()
        self.ymid = .5*(self.ystart+self.yend).flatten()
        self.zmid = .5*(self.zstart+self.zend).flatten()

    def set_imem_recorders(self):
        '''
        Record membrane currents for all segments
        '''
        self.memireclist = neuron.h.List()
        for sec in self.allseclist:
            for seg in sec:
                memirec = h.Vector(int(h.tstop/h.dt +1))
                memirec.record(seg._ref_i_membrane)
                self.memireclist.append(memirec)
    def set_voltage_recorders(self):
        '''
        Record voltage for all segments (not from LFPy sources)
        '''
        self.vreclist = h.List()

        for sec in self.allseclist:
            # address the problem of the important number of segments necessary to compute the accurate AP propagation in the unmyelinated axon case
            if sec.nseg > 100:
                for i in range(1,101):
                    vrec = h.Vector(int(h.tstop/h.dt+1))
                    vrec.record(sec(float(i)/100)._ref_v)
                    self.vreclist.append(vrec)
            else:
                for seg in sec:
                    vrec = h.Vector(int(h.tstop/h.dt+1))
                    vrec.record(seg._ref_v)
                    self.vreclist.append(vrec)


    def calc_imem(self):
        '''
        Fetch the vectors from the memireclist and calculate self.imem
        containing all the membrane currents.
        '''
        self.imem = np.array(self.memireclist)
        for i in range(self.imem.shape[0]):
            self.imem[i, ] *= self.area[i] * 1E-2
        self.memireclist = None
        del self.memireclist

    def collect_tvec(self):
        '''
        Set the tvec to be a monotonically increasing numpy array after sim.
        '''
        self.tvec = np.arange(h.tstop /h.dt + 1)*h.dt

    def simulate(self, electrode=None, rec_imem=True,
                 rec_variables=[], variable_dt=False, atol=0.001,
                 to_memory=True, to_file=False, file_name=None,
                 dotprodcoeffs=None):
        '''rec_vmem=False,
                 rec_ipas=False, rec_icap=False,
                 rec_isyn=False, rec_vmemsyn=False, rec_istim=False,'''
        '''
        This is the main function running the simulation of the NEURON model.
        Start NEURON simulation and record variables specified by arguments.

        Arguments:
        ::

            electrode:  Either an LFPy.RecExtElectrode object or a list of such.
                        If supplied, LFPs will be calculated at every time step
                        and accessible as electrode.LFP. If a list of objects
                        is given, accessible as electrode[0].LFP etc.
            rec_imem:   If true, segment membrane currents will be recorded
                        If no electrode argument is given, it is necessary to
                        set rec_imem=True in order to calculate LFP later on.
                        Units of (nA).
            rec_v:   record segment voltages (mV)
            rec_istim:  record currents of StimIntraElectrode (nA)
            rec_variables: list of variables to record, i.e arg=['cai', ]
            variable_dt: boolean, using variable timestep in NEURON
            atol:       absolute tolerance used with NEURON variable timestep
            to_memory:  only valid with electrode, store lfp in -> electrode.LFP
            to_file:    only valid with electrode, save LFPs in hdf5 file format
            file_name:  name of hdf5 file, '.h5' is appended if it doesnt exist
            dotprodcoeffs :  list of N x Nseg np.ndarray. These arrays will at
                        every timestep be multiplied by the membrane currents.
                        Presumably useful for memory efficient csd or lfp calcs
            '''
        self.collect_tvec()

        if rec_imem:
            self.set_imem_recorders()
        if self.rec_v:
            self.set_voltage_recorders()
 

        #run fadvance until t >= tstopms, and calculate LFP if asked for
        """if electrode is None and dotprodcoeffs is None:
            if not rec_imem:
                print(("rec_imem = %s, membrane currents will not be recorded!" \
                                  % str(rec_imem)))
            _run_simulation(self, variable_dt, atol)
        else:"""
        #allow using both electrode and additional coefficients:
        run_simulation._run_simulation_with_electrode(self, electrode, variable_dt, atol,
                                               to_memory, to_file, file_name,
                                               dotprodcoeffs)
        if rec_imem:
            self.calc_imem()

    """
    Equivalent methods interpxyz and setrx from the xtra mechanism available on the NEURON website from Ted Carnevale
    Setrx has been modified to integrate the use of multipolar electrodes
    """
    def interpxyz(self):
        """ interpolated data, spaced at regular intervals
        """
        # First, need to interpolate centers unto all compartments; from interpxyz.hoc
        for sec in self.allseclist:
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
                '''to use Vector class's .interpolate() must first scale the independent variable i.e. normalize length along centroid'''
                length.div(length.x[nn-1])
                # initialize the destination "independent" vector
                rr = h.Vector(sec.nseg+2)
                rr.indgen(1./sec.nseg)
                rr.sub(1./(2.*sec.nseg))
                rr.x[0]=0.
                rr.x[sec.nseg+1]=1.
                ''' length contains the normalized distances of the pt3d points along the centroid of the section.  
                These are spaced at irregular intervals. 
                range contains the normalized distances of the nodes along the centroid of the section.  
                These are spaced at regular intervals.'''
                # Ready to interpolate.
                xint = h.Vector(sec.nseg+2)
                yint = h.Vector(sec.nseg+2)
                zint = h.Vector(sec.nseg+2)
                xint.interpolate(rr, length, xx)
                yint.interpolate(rr, length, yy)
                zint.interpolate(rr, length, zz)
                '''for each node, assign the xyz values to x_xtra, y_xtra, z_xtra
                don't bother computing coords of the 0 and 1 ends
                also avoid writing coords of the 1 end into the last internal node's coords'''
                for ii in range(1,sec.nseg+1):
                    xr = rr.x[ii]
                    sec(xr).x_xtra = xint.x[ii]
                    sec(xr).y_xtra = yint.x[ii]
                    sec(xr).z_xtra = zint.x[ii]

    def setrx(self,stim_elec, axon_pos):
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
                        """ 0.01 converts rho's cm to um and ohm to megohm
                        if electrode is exactly at a node, r will be 0
                        this would be meaningless since the location would be inside the cell
                        so force r to be at least as big as local radius """
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
                        

class Unmyelinated(Axon):
    """
    name: axon name (for neuron)
    
    nsegs_method: ['lambda100']/'lambda_f'/'fixed_length': nseg rule
    max_nsegs_length: [None]: max segment length for method 'fixed_length'
    lambda_f: [100]: AC frequency for method 'lambda_f'
    d_lambda: [0.1]: parameter for d_lambda rule

    L: Axon length (micrometer)
    diam: Axon diameter (micrometer)
    cm : Specific membrane capacitance (microfarad/cm2)
    Ra: Specific axial resistance (Ohm cm)
    coord: y,z coordinates to spatially place the axon once created
    layout3D: either "DEFINE_SHAPE" or "PT3D" using hoc corresponding function
    rec_v: set voltage recorders True or False
    """
    def __init__(self, name, L, diam, cm, Ra, coord, layout3D, rec_v, hhDraw=False, nsegs_method='lambda100', lambda_f=100, d_lambda=0.1, max_nsegs_length=None):
        super(Unmyelinated,self).__init__(layout3D, rec_v)

        # set all properties
        self.coord = coord
        self.L = L
        self.diam = diam
        self.fiberD = diam
        print "Unmyelinated axon diameter: " +str(self.diam)
        self.cm = cm
        self.Ra = Ra
        self.hhDraw = hhDraw
        self.name = name
        self.nsegs_method = nsegs_method
        self.lambda_f = lambda_f
        self.d_lambda = d_lambda
        self.max_nsegs_length = max_nsegs_length
        ## End insertion

    def create_neuron_object(self):
        self.axon = h.Section(name = str(self.name))
        self.allseclist = h.SectionList()
        self.allseclist.append(sec = self.axon)

        self.axon.insert('extracellular')
        self.axon.insert('xtra')
        self.axon_update_property()
        if (self.layout3D == "DEFINE_SHAPE"):
            h.define_shape()
        elif (self.layout3D == "PT3D"):
             ### PLACE AXON SPATIALLY #######################################################        
            h.pt3dclear(sec=self.axon)
            h.pt3dadd(0, self.coord[0], self.coord[1], self.axon.diam, sec=self.axon)
            h.pt3dadd(self.axon.L, self.coord[0], self.coord[1], self.axon.diam, sec=self.axon)
            ################################################################################
        else:
            raise NameError('layout3D only "DEFINE_SHAPE" or "PT3D"')

        self.set_nsegs(self.nsegs_method, self.lambda_f, self.d_lambda, self.max_nsegs_length)

        for sec_id in self.allseclist:
            for seg in sec_id:
                h.setpointer(seg._ref_i_membrane, 'im', seg.xtra)
                h.setpointer(seg._ref_e_extracellular, 'ex', seg.xtra)
        self.interpxyz()
        self.collect_geometry()
        self.channel_init(0.120,0.036,0.0003,50,-77,-54.3) # default values of hh channel are used
        
    def delete_neuron_object(self):
        ### DELETE THE PYTHON REFERENCE TO THE AXON TO DELETE THE AXON ###
        for sec in self.allseclist:
            for seg in sec:
                seg = None
            sec = None
        self.axon = None

    def axon_update_property(self):
        self.axon.L = self.L
        self.axon.diam = self.diam
        self.axon.cm = self.cm
        self.axon.Ra = self.Ra
        
        
    """
    g unit S/cm2 and e unit mV
    gna: maximum specific sodium channel conductance
    gk: maximum specific potassium channel conductance
    gl: maximum specific leakage conductance
    ena: reversal potential for the sodium channel
    ek: reversal potential for the potassium channel
    el: reversal potential for the leakage channel
    """
    def channel_init(self, gna,gk,gl, ena,ek,el):
        self.axon.insert('hh') # insert a Hodgkin & Huxley channel
        if not self.hhDraw:
            self.axon.gnabar_hh = gna
            self.axon.gkbar_hh = gk
            self.axon.gl_hh = gl
            self.axon.ena = ena
            self.axon.ek = ek
            self.axon.el_hh = el
        else:
            self.axon.gnabar_hh = gna*(1+0.1*np.random.uniform(1))#[0])
            self.axon.gkbar_hh = gk*(1+0.1*np.random.uniform(1))#[0])
            self.axon.gl_hh = gl
            self.axon.ena = ena*(1+0.2*np.random.uniform(1))#[0])
            self.axon.ek = ek*(1+0.1*np.random.uniform(1))#[0])
            self.axon.el_hh = el

    ## LPFy methods using the d_lambda rule available in NEURON
    # http://www.neuron.yale.edu/neuron/static/docs/d_lambda/d_lambda.html

    def set_nsegs(self, nsegs_method, lambda_f, d_lambda, max_nsegs_length):
        '''Set number of segments per section according to the lambda-rule,
        or according to maximum length of segments'''
        if nsegs_method == 'lambda100':
            self.set_nsegs_lambda100(d_lambda)
        elif nsegs_method == 'lambda_f':
            self.set_nsegs_lambda_f(lambda_f, d_lambda)
        elif nsegs_method == 'fixed_length':
            self.set_nsegs_fixed_length(max_nsegs_length)
        else:
            if self.verbose:
                print(('No nsegs_method applied (%s)' % nsegs_method))
        self.totnsegs = self.calc_totnsegs()
                
    def set_nsegs_lambda_f(self, frequency=100, d_lambda=0.1):
        '''Set the number of segments for section according to the 
        d_lambda-rule for a given input frequency
            frequency: float, frequency at whihc AC length constant is computed
            d_lambda: float, 
        '''
        for sec in self.allseclist:
            sec.nseg = int((sec.L / (d_lambda*h.lambda_f(frequency,
                                                           sec=sec)) + .9)
                / 2 )*2 + 1
            print "Number of segments for unmyelinated axons via d_lambda: "+ str(sec.nseg)
        if self.verbose:
            print(("set nsegs using lambda-rule with frequency %i." % frequency))
   
    def set_nsegs_lambda100(self, d_lambda=0.1):
        '''Set the numbers of segments using d_lambda(100)'''
        self.set_nsegs_lambda_f(frequency=100, d_lambda=d_lambda)
    
    def set_nsegs_fixed_length(self, maxlength):
        '''Set nseg for sections so that every segment L < maxlength'''
        for sec in self.allseclist:
            sec.nseg = int(sec.L / maxlength) + 1







'''/*--------------------------------------------------------------------
2/02
Cameron C. McIntyre
SIMULATION OF PNS MYELINATED AXON

This model is described in detail in:

McIntyre CC, Richardson AG, and Grill WM. Modeling the excitability of
mammalian nerve fibers: influence of afterpotentials on the recovery
cycle. Journal of Neurophysiology 87:995-1006, 2002.
----------------------------------------------------------------------*/'''
class Myelinated(Axon):
    """
    name: axon name (for neuron)
    Nnodes: Number of nodes
    fiberD: Diameter of the fiber
    coord: y,z coordinates to spatially place the axon once created
    layout3D: either "DEFINE_SHAPE" or "PT3D" using hoc corresponding function
    rec_v: set voltage recorders True or False
    ### PARAMETERS BELOW ARE ACTUALLY SET IN THE CODE BUT USEFUL FOR THE HEADER ###
    ## Values of paralength2 and interlength are to be set according to the chosen fiberD value
    nodelength : set node length
    paralength1: set the length of the nearest paranode segment to the node
    paralength2: set the length of the second paranode segment followed by the internodes segments
    interlength: set the length of the internode part comprising the 6 segments between two paranodes2
    """
    def __init__(self, name, Nnodes, fiberD, coord, layout3D, rec_v):#, nodelength, paralength1, paralength2, interlength):
        
        self.name = name
        self.layout3D = layout3D
        self.axonnodes = Nnodes
        self.coord = coord
        self.fiberD = fiberD #choose from 5.7, 7.3, 8.7, 10.0, 11.5, 12.8, 14.0, 15.0, 16.0

        print 'Myelinated fiber diameter: ' + str(self.fiberD)
        # topological parameters
        """self.axonnodes=21"""
        self.paranodes1= 2*(self.axonnodes-1)
        self.paranodes2= 2*(self.axonnodes-1)
        self.axoninter= 6*(self.axonnodes-1)
        self.axontotal= self.axonnodes+self.paranodes1+self.paranodes2+self.axoninter

        # morphological parameters
        self.paralength1=3
        self.nodelength=1.0
        self.space_p1=0.002
        self.space_p2=0.004
        self.space_i=0.004

        #electrical parameters
        self.rhoa=0.7e6 #Ohm-um
        self.mycm=0.1 #uF/cm2/lamella membrane
        self.mygm=0.001 #S/cm2/lamella membrane



        ###### from proc dependent_var() #############
        if (self.fiberD==5.7):
            g=0.605
            axonD=3.4
            nodeD=1.9
            paraD1=1.9
            paraD2=3.4
            deltax=500
            self.paralength2=35
            nl=80
        if (self.fiberD==7.3):
            g=0.630
            axonD=4.6
            nodeD=2.4
            paraD1=2.4
            paraD2=4.6
            deltax=750
            self.paralength2=38
            nl=100
        if (self.fiberD==8.7):
            g=0.661
            axonD=5.8
            nodeD=2.8
            paraD1=2.8
            paraD2=5.8
            deltax=1000
            self.paralength2=40
            nl=110
        if (self.fiberD==10.0):
            g=0.690
            axonD=6.9
            nodeD=3.3
            paraD1=3.3
            paraD2=6.9
            deltax=1150
            self.paralength2=46
            nl=120
        if (self.fiberD==11.5):
            g=0.700
            axonD=8.1
            nodeD=3.7
            paraD1=3.7
            paraD2=8.1
            deltax=1250
            self.paralength2=50
            nl=130
        if (self.fiberD==12.8):
            g=0.719
            axonD=9.2
            nodeD=4.2
            paraD1=4.2
            paraD2=9.2
            deltax=1350
            self.paralength2=54
            nl=135
        if (self.fiberD==14.0):
            g=0.739
            axonD=10.4
            nodeD=4.7
            paraD1=4.7
            paraD2=10.4
            deltax=1400
            self.paralength2=56
            nl=140
        if (self.fiberD==15.0):
            g=0.767
            axonD=11.5
            nodeD=5.0
            paraD1=5.0
            paraD2=11.5
            deltax=1450
            self.paralength2=58
            nl=145
        if (self.fiberD==16.0):
            g=0.791
            axonD=12.7
            nodeD=5.5
            paraD1=5.5
            paraD2=12.7
            deltax=1500
            self.paralength2=60
            nl=150

        self.g = g
        self.axonD = axonD
        self.nodeD = nodeD
        self.paraD1 = paraD1
        self.paraD2 = paraD2
        self.deltax = deltax
        self.nl = nl

        self.Rpn0=(self.rhoa*0.01)/(math.pi*((math.pow((self.nodeD/2)+self.space_p1,2))-(math.pow(nodeD/2,2))))
        self.Rpn1=(self.rhoa*0.01)/(math.pi*((math.pow((self.paraD1/2)+self.space_p1,2))-(math.pow(paraD1/2,2))))
        self.Rpn2=(self.rhoa*0.01)/(math.pi*((math.pow((self.paraD2/2)+self.space_p2,2))-(math.pow(paraD2/2,2))))
        self.Rpx=(self.rhoa*0.01)/(math.pi*((math.pow((self.axonD/2)+self.space_i,2))-(math.pow(axonD/2,2))))

        self.interlength=(self.deltax-self.nodelength-(2*self.paralength1)-(2*self.paralength2))/6

        super(Myelinated,self).__init__(layout3D, rec_v)

    def create_neuron_object(self):

        self.allseclist = h.SectionList()

        #### from initialize() ####
        self.nodes = []
        for i in range(self.axonnodes):
            node = h.Section()
            node.nseg = 1
            node.diam = self.nodeD
            node.L = self.nodelength
            node.Ra = self.rhoa/10000
            node.cm = 2
            node.insert('axnode')
            #h('insert extracellular xraxial=Rpn0 xg=1e10 xc=0')
            node.insert('extracellular')
            node.xraxial[0] = self.Rpn0
            node.xraxial[1] = self.Rpn0
            node.xg[0] =1e10
            node.xg[1] =1e10
            node.xc[0] =0
            node.xc[1] =0
            node.insert('xtra')

            self.nodes.append(node)
            self.allseclist.append(sec=node)

        self.MYSAs = [] # paranodes1
        for i in range(self.paranodes1):
            MYSA = h.Section()
            MYSA.nseg = 1
            MYSA.diam = self.fiberD
            MYSA.L = self.paralength1
            MYSA.Ra = self.rhoa*(1/math.pow(self.paraD1/self.fiberD,2))/10000
            MYSA.cm = 2*self.paraD1/self.fiberD
            MYSA.insert('pas')
            MYSA.g_pas = 0.001*self.paraD1/self.fiberD
            MYSA.e_pas = -80
            #h('insert extracellular xraxial=Rpn1 xg=mygm/(nl*2) xc=mycm/(nl*2)')
            MYSA.insert('extracellular')
            MYSA.xraxial[0] = self.Rpn1
            MYSA.xraxial[1] = self.Rpn1
            MYSA.xg[0] = self.mygm/(self.nl*2)
            MYSA.xg[1] = self.mygm/(self.nl*2)
            MYSA.xc[0] = self.mycm/(self.nl*2)
            MYSA.xc[1] = self.mycm/(self.nl*2)
            MYSA.insert('xtra')

            self.MYSAs.append(MYSA)
            self.allseclist.append(sec=MYSA)

        self.FLUTs = [] # paranodes2
        for i in range(self.paranodes2):
            FLUT = h.Section()
            FLUT.nseg = 1
            FLUT.diam = self.fiberD
            FLUT.L = self.paralength2
            FLUT.Ra = self.rhoa*(1/math.pow(self.paraD2/self.fiberD,2))/10000
            FLUT.cm = 2*self.paraD2/self.fiberD
            FLUT.insert('pas')
            FLUT.g_pas = 0.0001*self.paraD2/self.fiberD
            FLUT.e_pas = -80
            #h('insert extracellular xraxial=Rpn2 xg=mygm/(nl*2) xc=mycm/(nl*2)')
            FLUT.insert('extracellular')
            FLUT.xraxial[0] = self.Rpn2
            FLUT.xraxial[1] = self.Rpn2
            FLUT.xg[0] = self.mygm/(self.nl*2)
            FLUT.xg[1] = self.mygm/(self.nl*2)
            FLUT.xc[0] = self.mycm/(self.nl*2)
            FLUT.xc[1] = self.mycm/(self.nl*2)
            FLUT.insert('xtra')

            self.FLUTs.append(FLUT)
            self.allseclist.append(sec=FLUT)

        self.STINs = [] # internodes
        for i in range(self.axoninter):
            STIN = h.Section()
            STIN.nseg = 1
            STIN.diam = self.fiberD
            STIN.L = self.interlength
            STIN.Ra = self.rhoa*(1/math.pow(self.axonD/self.fiberD,2))/10000
            STIN.cm = 2*self.axonD/self.fiberD
            STIN.insert('pas')
            STIN.g_pas = 0.0001*self.axonD/self.fiberD
            STIN.e_pas = -80
            #h('insert extracellular xraxial=Rpx xg=mygm/(nl*2) xc=mycm/(nl*2)')
            STIN.insert('extracellular')
            STIN.xraxial[0] = self.Rpx
            STIN.xraxial[1] = self.Rpx
            STIN.xg[0] = self.mygm/(self.nl*2)
            STIN.xg[1] = self.mygm/(self.nl*2)
            STIN.xc[0] = self.mycm/(self.nl*2)
            STIN.xc[1] = self.mycm/(self.nl*2)
            STIN.insert('xtra')

            self.STINs.append(STIN)
            self.allseclist.append(sec=STIN)

        if (self.layout3D == "DEFINE_SHAPE"):
            for i in range(self.axonnodes-1):
                self.MYSAs[2*i].connect(self.nodes[i],1,0)
                self.FLUTs[2*i].connect(self.MYSAs[2*i],1,0)
                self.STINs[6*i].connect(self.FLUTs[2*i],1,0)
                self.STINs[6*i+1].connect(self.STINs[6*i],1,0)
                self.STINs[6*i+2].connect(self.STINs[6*i+1],1,0)
                self.STINs[6*i+3].connect(self.STINs[6*i+2],1,0)
                self.STINs[6*i+4].connect(self.STINs[6*i+3],1,0)
                self.STINs[6*i+5].connect(self.STINs[6*i+4],1,0)
                self.FLUTs[2*i+1].connect(self.STINs[6*i+5],1,0)
                self.MYSAs[2*i+1].connect(self.FLUTs[2*i+1],1,0)
                self.nodes[i+1].connect(self.MYSAs[2*i+1],1,0)
            h.define_shape()
        elif (self.layout3D == "PT3D"):
            ### PLACE AXON SPATIALLY #############################################################################################################################################################       
            for i in range(self.axonnodes-1):
                h.pt3dclear(sec=self.nodes[i])
                h.pt3dadd(0, self.coord[0], self.coord[1], self.nodes[i].diam, sec=self.nodes[i])
                h.pt3dadd(self.nodelength*i+self.interlength*6*i+2*i*self.paralength1+ 2*i*self.paralength2 + self.nodelength, self.coord[0], self.coord[1], self.nodes[i].diam, sec=self.nodes[i])

                self.MYSAs[2*i].connect(self.nodes[i],1,0)
                h.pt3dclear(sec=self.MYSAs[2*i])
                h.pt3dadd(self.nodelength*i+self.interlength*6*i+2*i*self.paralength1+ 2*i*self.paralength2 + self.nodelength, self.coord[0], self.coord[1], self.MYSAs[2*i].diam, sec=self.MYSAs[2*i])
                h.pt3dadd(self.nodelength*(i+1)+self.interlength*6*i+2*i*self.paralength1+  2*i*self.paralength2 +self.paralength1, self.coord[0], self.coord[1], self.MYSAs[2*i].diam, sec=self.MYSAs[2*i])

                self.FLUTs[2*i].connect(self.MYSAs[2*i],1,0)
                h.pt3dclear(sec=self.FLUTs[2*i])
                h.pt3dadd(self.nodelength*(i+1)+self.interlength*6*i+2*i*self.paralength1+  2*i*self.paralength2 +self.paralength1, self.coord[0], self.coord[1], self.FLUTs[2*i].diam, sec=self.FLUTs[2*i])
                h.pt3dadd(self.nodelength*(i+1)+self.interlength*6*i+2*i*self.paralength1+  2*i*self.paralength2 + self.paralength1 + self.paralength2, self.coord[0], self.coord[1], self.FLUTs[2*i].diam, sec=self.FLUTs[2*i])
                for j in range(6):
                    self.STINs[6*i+j].connect(self.FLUTs[2*i],1,0)
                    h.pt3dclear(sec=self.STINs[6*i+j])
                    h.pt3dadd(self.nodelength*(i+1)+self.interlength*6*i+2*i*self.paralength1+  2*i*self.paralength2 + self.paralength1 + self.paralength2+ j*self.interlength, self.coord[0], self.coord[1], self.STINs[6*i+j].diam, sec=self.STINs[6*i+j])
                    h.pt3dadd(self.nodelength*(i+1)+self.interlength*6*i+2*(i+1.0/2)*self.paralength1+  2*(i+1.0/2)*self.paralength2 + (j+1)*self.interlength, self.coord[0], self.coord[1], self.STINs[6*i+j].diam, sec=self.STINs[6*i+j])
                self.FLUTs[2*i+1].connect(self.STINs[6*i+5],1,0)
                h.pt3dclear(sec=self.FLUTs[2*i+1])
                h.pt3dadd(self.nodelength*(i+1)+self.interlength*6*i+2*(i+1.0/2)*self.paralength1+  2*(i+1.0/2)*self.paralength2 + 6*self.interlength, self.coord[0], self.coord[1], self.FLUTs[2*i+1].diam, sec=self.FLUTs[2*i+1])
                h.pt3dadd(self.nodelength*(i+1)+self.interlength*6*(i+1)+2*(i+1.0/2)*self.paralength1+  2*(i+1.0/2)*self.paralength2 +self.paralength1, self.coord[0], self.coord[1], self.FLUTs[2*i+1].diam, sec=self.FLUTs[2*i+1])
                self.MYSAs[2*i+1].connect(self.FLUTs[2*i+1],1,0)
                h.pt3dclear(sec=self.MYSAs[2*i+1])
                h.pt3dadd(self.nodelength*(i+1)+self.interlength*6*(i+1)+2*(i+1.0/2)*self.paralength1+  2*(i+1.0/2)*self.paralength2 +self.paralength1, self.coord[0], self.coord[1], self.MYSAs[2*i+1].diam, sec=self.MYSAs[2*i+1])
                h.pt3dadd(self.nodelength*(i+1)+self.interlength*6*(i+1)+2*(i+1.0/2)*self.paralength1+  2*(i+1.0/2)*self.paralength2 + self.paralength1 + self.paralength2, self.coord[0], self.coord[1], self.MYSAs[2*i+1].diam, sec=self.MYSAs[2*i+1])
                self.nodes[i+1].connect(self.MYSAs[2*i+1],1,0)

            h.pt3dclear(sec=self.nodes[self.axonnodes-1])
            h.pt3dadd(self.nodelength*(self.axonnodes-1)+self.interlength*(self.axonnodes-1)*6+2*(self.axonnodes-1)*self.paralength1+ 2*(self.axonnodes-1)*self.paralength2, self.coord[0], self.coord[1], self.nodes[self.axonnodes-1].diam, sec=self.nodes[self.axonnodes-1])
            h.pt3dadd(self.nodelength*(self.axonnodes-1)+self.interlength*(self.axonnodes-1)*6+2*(self.axonnodes-1)*self.paralength1+ 2*(self.axonnodes-1)*self.paralength2 +self.nodelength, self.coord[0], self.coord[1], self.nodes[self.axonnodes-1].diam, sec=self.nodes[self.axonnodes-1])
        #######################################################################################################################################################################################################
        else:
            raise NameError('layout3D only "DEFINE_SHAPE" or "PT3D"')


        self.totnsegs = self.calc_totnsegs()
        for sec_id in self.allseclist:
            for seg in sec_id:
                h.setpointer(seg._ref_i_membrane, 'im', seg.xtra)
                h.setpointer(seg._ref_e_extracellular, 'ex', seg.xtra)
        self.interpxyz()
        self.collect_geometry()
        
    def delete_neuron_object(self):
        ### DELETE THE PYTHON REFERENCE TO THE AXON TO DELETE THE AXON ###
        for sec in self.allseclist:
            for seg in sec:
                seg = None
            sec = None

        self.nodes = None
        self.FLUTs = None
        self.MYSAs = None
        self.STINs = None

class Stimulus(object):
    """
    stim_type: INTRA or EXTRA cellular stimulation
    axon: axon object on which the stimulation is applied
    pos: position of the stimulus
    sect: section being stimulated
    delay: pulse delay (ms)
    dur: pulse duration (ms)
    amp: pulse amplitude (nA)
    freq: frequency of the sin pulse (Hz)
    duty_cycle: Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
    stim_coord=[xe,ye,ze]: spatial coordinates  of the stimulating electrode
    waveform: Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
    """
    def __init__(self, stim_type, dur,amp, freq, duty_cycle, stim_coord, waveform):
        self.waveform = waveform
        self.stim_type = stim_type
        self.stim_coord = stim_coord
        self.dur = dur
        self.freq = freq
        self.amp = amp
        self.duty_cycle = duty_cycle
        cut_off = math.sin((-self.duty_cycle+1.0/2)*math.pi)
        self.t = np.linspace(0, self.dur, (h.tstop+2*h.dt)/h.dt, endpoint=True)
        if self.waveform == "MONOPHASIC":
            self.signal = self.amp*(cut_off < np.sin(2*math.pi*self.freq*(self.t)))
        elif self.waveform == "BIPHASIC":
            self.signal = self.amp*(cut_off < np.sin(2*math.pi*self.freq*self.t))-self.amp*(cut_off < np.sin(2*math.pi*self.freq*(self.t+self.duty_cycle/self.freq)))
        else:
            print "You didn't choose the right waveform either MONOPHASIC or BIPHASIC, it has been set to default MONOPHASIC"
            self.signal = self.amp*(cut_off < np.sin(2*math.pi*self.freq*(self.t)))
            
        self.svec = h.Vector(self.signal)

    def connectAxon(self, axon):
        if self.stim_type == "INTRA":
            self.init_intra(axon)
        elif self.stim_type == "EXTRA":
            axon.setrx(self.stim_coord, axon.coord)
            self.init_xtra()
        else:
            raise NameError('stim_type only "INTRA" or "EXTRA"')

    def init_intra(self, axon):
        # Place an IClamp on the first element of the allseclist
        # In unmyelinated axon case allseclist is directly the unique axon section

        # counter = 0
        # for seg in axon.allseclist[0]:
        #     print 'Segment number ' + str(counter)
        #     print seg
        #     counter += 1

        axon.stim = h.IClamp(0, axon.allseclist)
        axon.stim.delay = 0
        axon.stim.dur = self.dur
        self.svec.play(axon.stim._ref_amp,h.dt)
        
    def init_xtra(self):
        self.svec.play(h._ref_is_xtra,h.dt)

        
class Bundle(object):
    """
    radius_bundle: Radius of the bundle (typically 0.5-1.5mm)
    number_of_axons: Number of axons in the bundle
    p_A: Percentage of myelinated fiber type A
    p_B: Percentage of myelinated fiber type B
    p_C: Percentage of unmyelinated fiber type C
    number_contact_points: Number of points on the circle constituing the cuff electrode
    recording_elec: Position of the recording electrode along axon in um, in "BIPOLAR" case the position along axons should be given as a couple [x1,x2]
    jitter_para: Mean and standard deviation of the jitter
    stim_type: Stimulation type either "INTRA" or "EXTRA" for INTRA/EXTRA_cellular stimulation
    stim_coord: spatial coordinates  of the stimulating electrode, example for tripolar case=[[xe0,ye0,ze0], [xe1,ye1,ze1], [xe2,ye2,ze2]] (order is important with the middle being the cathode)
    amplitude: pulse amplitude (nA)
    freq: frequency of the sin pulse (kHz)
    duty_cycle: Percentage stimulus is ON for one cycl
    stim_dur : stimulus duration (ms)
    dur: simulation duration (ms)
    waveform: Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
    number_of_elecs: number of electrodes along the bundle
    recording_type: "BIPOLAR" or "MONOPOLAR"
    myelinated_A: parameters for fiber type A
    umyelinated:  parameters for fiber type C
    draw_distribution: Boolean stating the distribution of fibre should be drawn

    
    """
    def __init__(self, radius_bundle, draw_distribution, number_of_axons, p_A, p_C, number_contact_points, recording_elec_pos, jitter_para, stim_type, stim_coord, duty_cycle, freq, amplitude, stim_dur, dur, number_elecs, myelinated_A, unmyelinated,rec_CAP, waveform):

        self.draw_distribution = draw_distribution
        self.myelinated_A =  myelinated_A
        self.unmyelinated =  unmyelinated

        self.waveform = waveform
        self.stim_type = stim_type
        #self.stim_coord = stim_coord#[0,500,0] # [xe,ye,ze] #um
        self.stim_coord = [[0, radius_bundle*math.cos(math.pi/number_contact_points), radius_bundle*math.sin(math.pi/number_contact_points)]]
        self.duty_cycle = duty_cycle
        self.freq = freq
        self.amp = amplitude
        self.stim_dur = stim_dur
        self.dur = dur
        self.rec_CAP = rec_CAP

        self.p_A = float(p_A)/(p_A+p_C)
        self.p_C = float(p_C)/(p_A+p_C)

        self.number_of_axons = number_of_axons
        self.axons = []
        self.radius_bundle = radius_bundle # um
        self.electrodes = []
        self.voltages = []
        self.number_contact_points = number_contact_points
        self.number_elecs = number_elecs
        self.recording_elec_pos = recording_elec_pos #um
       
        self.build_disk(self.number_of_axons,self.radius_bundle)

        self.saveParams={'elecCount': len(self.recording_elec_pos), 'dt': h.dt, 'tStop': h.tstop, 'p_A': self.p_A,
                    'myelinatedDiam': self.myelinated_A['fiberD'], 'unmyelinatedDiam': self.unmyelinated['diam'],
                    'L': self.unmyelinated['L'], 'stimType': self.stim_type, 'stimWaveform' : self.waveform,
                    'stimDutyCycle': self.duty_cycle, 'stimAmplitude' : self.amp}

        # ### JITTER (random gaussian delay for individual fiber stimulation) ###
        # self.delay_mu, self.delay_sigma = jitter_para[0], jitter_para[1] # mean and standard deviation
        # if (jitter_para[0] == 0) & (jitter_para[1] == 0):
        #     delay = np.zeros(self.number_of_axons)
        # else:
        #     delay = np.random.normal(self.delay_mu, self.delay_sigma, self.number_of_axons)

        # create axons within wanted nerve-slice


        # create axons within one fourth slice of the whole bundle.
        self.virtual_number_axons = 0
        for i in range(self.number_of_axons):
            # if within the desired subsection of the bundle, create axon
            if ((self.axons_pos[i,1]>=0 and self.axons_pos[i,1]< self.axons_pos[i,0]) or (self.axons_pos[i,0]== 0 and self.axons_pos[i,1]==0)):
                # print self.axons_pos[i,:]
                self.create_axon(self.axons_pos[i,:])
                self.virtual_number_axons +=1
                print "Number axons created:" + str(self.virtual_number_axons)

        # create Simulus instace used for all axons
        self.stim = Stimulus(self.stim_type, self.stim_dur,self.amp, self.freq,self.duty_cycle, self.stim_coord, self.waveform)

        # # connect to stimulus
        # for i in range(self.virtual_number_axons):
        #     axon = self.axons[i]
        #     self.stim.connectAxon(axon)

    def simulateBundle(self):

        self.simulateAxons()

        if self.rec_CAP:
            if self.number_elecs !=1:
                self.compute_CAP2D_fromfile() 
            else:
                self.compute_CAP1D_fromfile()
            self.save_CAP_to_file()

        self.save_voltage_to_file()
                    
    def save_CAP_to_file(self):
        filename = self.get_filename("CAP")
        print "Save location for CAP file: " + filename

        # maybe add the header later. Right now we assume that the simulation is defined by the bundle object that get
        # always generated during the whole simulation. If files are opened independently from a bundle object, such a
        # header would be useful.
        # header = repr(parameters)
        DataOut = np.array(self.trec)
        if self.number_elecs != 1:
            for i in range(len(self.sum_CAP)):
                DataOut = np.column_stack( (DataOut, np.array(self.sum_CAP[i])))
        else:
            DataOut = np.column_stack( (DataOut, np.array(self.sum_CAP)))

        np.savetxt(filename, DataOut)

    def save_voltage_to_file(self):
        filename = self.get_filename("V")
        print "Save location for voltage file: " + filename
        #header= repr(parameters)

        DataOut = np.concatenate(([0],np.array(self.trec)))

        voltages = np.array(self.voltages)

        if np.size(voltages) == 0:
            return

        # as np.array(...) only converts outer Vector to python-readable format, we need to iterate through elements to convert
        # inner vectors where the actual voltage signals are stored.

        for i in range(len(voltages)):
            voltageSingleAxon = np.transpose(np.array(voltages[i]))

            # append the sectionlength in the first column in order to differentiate different axons later
            numberOfSegments = np.shape(voltageSingleAxon)[1]
            numberOfSegmentsArray = np.multiply(np.ones(numberOfSegments),np.array(numberOfSegments))
            voltageSingleAxonFormatted = np.row_stack((numberOfSegmentsArray, voltageSingleAxon))

            DataOut = np.column_stack( (DataOut, voltageSingleAxonFormatted))
        np.savetxt(filename, DataOut)#, header=header)

    def get_CAP_from_file(self):

        # get the whole CAP, can be single electrode or multiple
        directory = getDirectoryName("CAP", **self.saveParams)
        try:
            newestFile = max(glob.iglob(directory+'*.[Dd][Aa][Tt]'), key=os.path.getctime)
        except ValueError:
            print 'No CAP calculation has been performed yet with this set of parameter.'
            quit()

        CAPraw = np.transpose(np.loadtxt(newestFile))
        time = CAPraw[0,:]
        CAP = CAPraw[1:,:]

        return time, CAP

    def plot_CAP1D(self, maxNumberOfSubplots = 10):

        # first load the desired data from file
        time, CAP = self.get_CAP_from_file()

        numberOfRecordingSites = np.shape(CAP)[0]

        if not numberOfRecordingSites == 1:

            numberOfPlots = min(maxNumberOfSubplots, numberOfRecordingSites-1)

            eletrodeSelection = np.floor(np.linspace(1,numberOfRecordingSites-1, numberOfPlots))

            # Subplots
            f, axarr = plt.subplots(numberOfPlots, sharex=True)

            for i in range(numberOfPlots):

                electrodeIndex = eletrodeSelection[i]

                CAPSingleElectrode =  CAP[electrodeIndex,:]
                distanceFromOrigin = self.saveParams['L']/numberOfRecordingSites*electrodeIndex

                axarr[i].plot(time, CAPSingleElectrode)
                axarr[i].set_title('distance ' + str(distanceFromOrigin) + ' [um]')
                axarr[i].set_ylabel('CAP [mV]')

                if i == numberOfPlots - 1:
                    axarr[i].set_xlabel('time [ms]')

    def plot_CAP2D(self):

        # first load the desired data from file
        time, CAP = self.get_CAP_from_file()

        # print as an image
        fig = plt.figure()
        im = plt.imshow(CAP, cmap=plt.get_cmap('gist_stern'), interpolation='none', aspect='auto')#, norm=LogNorm(vmin=CAPmin, vmax=CAPmax))

        # correct xticks (from samples to ms)
        numberOfXTicks = 10
        tick_locs = np.round(np.linspace(0,np.shape(CAP)[1],numberOfXTicks))
        tick_lbls = np.round(np.linspace(0,self.saveParams['tStop'],numberOfXTicks))
        plt.xticks(tick_locs, tick_lbls, fontsize=12)

        # correct yticks (from electrodes to distance)
        numberOfYTicks = 10
        tick_locs = np.round(np.linspace(0,np.shape(CAP)[0],numberOfYTicks))
        tick_lbls = np.round(np.linspace(0,self.saveParams['L'],numberOfYTicks))
        plt.yticks(tick_locs, tick_lbls, fontsize=12)

        # add titles, axis labels and colorbar
        fig.suptitle('Compound action potential [uV] over space and time', fontsize=20)
        plt.xlabel('time [ms]')
        plt.ylabel('Distance from axon origin [um]')
        cbar = plt.colorbar(im)

    def plot_voltage(self):
        # get the whole CAP, can be signle electrode or multiple
        directory = getDirectoryName("V", **self.saveParams)
        try:
            newestFile = max(glob.iglob(directory+'*.[Dd][Aa][Tt]'), key=os.path.getctime)
        except ValueError:
            print 'No voltage calculation has been performed yet with this set of parameter.'
            quit()

        # load the raw voltage file
        timeStart = time.time()
        Vraw = np.transpose(np.loadtxt(newestFile))
        print 'Elapsed time to load voltage file ' + str(time.time() - timeStart) + 's'

        timeRec = Vraw[0,1:] # extract time vector
        segmentArray = Vraw[1:,0] # extract segment numbers for each axon (varies with diameter, lambda rule)
        V = Vraw[1:,1:] # free actual voltage signals from surrounding formatting

        # separate the axons, first get segment counts for each axon
        segmentNumbers = [segmentArray[0]]
        indexNextFirstEntry = segmentNumbers[0]
        while indexNextFirstEntry < len(segmentArray):
            segmentNumbers.append(segmentArray[indexNextFirstEntry])
            indexNextFirstEntry += segmentNumbers[-1]

        numberOfAxons = len(segmentNumbers)

        firstIndices = np.cumsum((np.insert(segmentNumbers,0,0)))

        voltageMatrices = []
        for i in range(numberOfAxons):
            startIndex = firstIndices[i]
            endIndex = firstIndices[i+1]-1
            voltageMatrices.append(V[startIndex:endIndex,:])

        # now plot
        numberOfAxons = np.shape(voltageMatrices)[0]
        numberOfPlots = min(6, numberOfAxons)

        axonSelection = np.floor(np.linspace(0,numberOfAxons-1, numberOfPlots))

        f, axarr = plt.subplots(numberOfPlots, sharex=True)

        # colors
        jet = plt.get_cmap('jet')

        for i in range(len(axonSelection)):
            axonIndex = int(axonSelection[i])

            voltageMatrix = np.transpose(voltageMatrices[axonIndex])

            # find out whether axon is myelinated or not
            isMyelinated = (type(self.axons[axonIndex]) == Myelinated)

            axonDiameter = self.axons[axonIndex].fiberD
            currentNumberOfSegments = np.shape(voltageMatrix)[1]

            if not isMyelinated:
                cNorm = colors.Normalize(vmin=0, vmax=currentNumberOfSegments-1)#len(diameters_m)-1)#
                scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)
                for j in range(currentNumberOfSegments):
                    colorVal = scalarMap.to_rgba(j)
                    axarr[i].plot(timeRec, voltageMatrix[:,j], color=colorVal)
                # axarr[i].set_title('distance ' + str(555) + ' [um]')
                axarr[i].set_ylabel('Voltage [mV]')
                axarr[i].set_xlabel('time [ms]')
                axarr[i].set_title('Voltage of unmyelinated axon with diameter ' + str(axonDiameter) + ' um')
            else:
                Nnodes = self.myelinated_A['Nnodes']

                cNorm = colors.Normalize(vmin=0, vmax=Nnodes-1)#len(diameters_m)-1)#
                scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

                for j in range(Nnodes):
                    colorVal = scalarMap.to_rgba(j)
                    axarr[i].plot(np.array(timeRec), np.array(voltageMatrix[:,j]), color=colorVal)

                axarr[i].set_ylabel('Voltage [mV]')
                axarr[i].set_xlabel('time [ms]')
                axarr[i].set_title('Voltage of nodes of myelinated axon with diameter ' + str(axonDiameter) + ' um')

    def create_axon(self, axonPosition):

        # first decide by chance whether myelinated or unmyelinated
        axonTypeIndex = np.random.choice(2,1,p = [self.p_A, self.p_C])
        axonTypes = ['m', 'u']
        axonType = axonTypes[axonTypeIndex]

        # then get diameter. Either drawn from distribution or constant.
        axonDiameter = self.getDiam(axonType)

        if axonTypeIndex == 1:
            unmyel = copy.copy(self.unmyelinated)
            unmyel['diam'] = axonDiameter
            axonParameters = dict( {'coord': axonPosition},**unmyel)
            self.axons.append(Unmyelinated(**axonParameters))

            # print 'before appending axon'
            # for i in range(len(self.axons)):
            #     print list(list(self.axons[i].allseclist)[0].allseg())
            #
            # newAxon = Unmyelinated(**axonParameters )
            #
            # print 'Segments of new axon:'
            # print list(list(newAxon.allseclist)[0].allseg())
            #
            # self.axons.append(newAxon)
            #
            # print 'after appending axon'
            # for i in range(len(self.axons)):
            #     print list(list(self.axons[i].allseclist)[0].allseg())

        elif axonTypeIndex == 0:
            myel = copy.copy(self.myelinated_A)
            myel['fiberD'] = axonDiameter
            axonParameters = dict( {'coord':axonPosition},**myel)
            self.axons.append(Myelinated(**axonParameters))
        else:
            "Error in the draw of the axon type!"

        # self.stim = Stimulus(self.stim_type,self.axons[i], delay[i],self.stim_dur,self.amp, self.freq,self.duty_cycle, self.stim_coord, self.waveform)

    def simulateAxons(self):

        # where are the electrodes
        [angles,X,Y,Z,N] = self.setup_recording_elec()

        runFlag = False
        for axonIndex in range(self.virtual_number_axons):

            axon = self.axons[axonIndex]

            # where is the axon
            axonPosition = axon.coord

            electrodeParameters = {         #parameters for RecExtElectrode class
                    'sigma' : 0.3,              #Extracellular potential
                    'x' : X,  #Coordinates of electrode contacts
                    'y' : Y-axonPosition[0],
                    'z' : Z-axonPosition[1],
                    'n' : 20,
                    'r' : 10,
                    'N' : N,
                    'method': "pointsource", #or "linesource"
                }

            ### LAUNCH SIMULATION ###
            temp = time.time()

            # create the neuron object specified in the axon class object
            axon.create_neuron_object()

            # connect to stimulus
            self.stim.connectAxon(axon)

            if axonIndex == 0:
            # record time variable
                self.trec = h.Vector()
                self.trec.record(h._ref_t)


            if self.rec_CAP:
                axon.simulate()
                self.electrodes.append(refextelectrode.RecExtElectrode(axon, **electrodeParameters))
                elapsed1 = time.time()-temp
                print "Elapsed time to simulate CAP: " + str(elapsed1)
                temp = time.time()
                self.electrodes[axonIndex].calc_lfp()
                elapsed2 = time.time()-temp
                print "Elapsed time to calc lfp:" + str(elapsed2)
                self.save_electrode(axonIndex)
                self.electrodes[axonIndex]= None
                self.CAP_to_file = True

            elif axon.rec_v:
                # just run a normal NEURON simulation to record the voltage
                axon.set_voltage_recorders()
                h.run() #runFlag = True #
                self.voltages.append(axon.vreclist)
                elapsed = time.time()-temp
                print "Elapsed time to simulate V: " + str(elapsed)
            else:
                print "You are probably recording nothing for this axon"
                h.run() #runFlag = True #

            # delete the object
            axon.delete_neuron_object()

        # if runFlag:
        #     startTime = time.time()
        #     h.run()
        #     print 'Elapsed time for voltage simulation ' + str(time.time() - startTime)
            # ### DELETE THE PYTHON REFERENCE TO THE AXON TO DELETE THE AXON ###
            # for sec in h.allsec():
            #     for seg in sec:
            #         seg = None
            #     sec = None
            # if draw[i] == 1:
            #     self.axons[i].axon = None
            # else:
            #     self.axons[i].nodes = None
            #     self.axons[i].FLUTs = None
            #     self.axons[i].MYSAs = None
            #     self.axons[i].STINs = None


    def store_geometry(self):
        self.geometry_parameters = [self.axons[0].xstart,self.axons[0].ystart,self.axons[0].zstart,self.axons[0].xend,self.axons[0].yend,self.axons[0].zend,self.axons[0].area,self.axons[0].diam,self.axons[0].length,self.axons[0].xmid,self.axons[0].ymid,self.axons[0].zmid]
    def save_electrode(self,i):
        directory = getDirectoryName("elec", **self.saveParams)

        print "saving electrode: "+str(i)
        if i==0:
            if not os.path.exists(directory):
                os.makedirs(directory)
            else:
                shutil.rmtree(directory)
                os.makedirs(directory)
        filename = "electrode_"+str(i)+".dat"
        DataOut = np.array(self.electrodes[i].LFP[0])
        for j in range(len(self.electrodes[i].LFP)):
            DataOut = np.column_stack((DataOut, np.array(self.electrodes[i].LFP[j])))
        np.savetxt(directory+filename, DataOut)
        
    def load_electrodes(self):
        directory = getDirectoryName("elec", **self.saveParams)

        print "loading electrode"
        t0 = time.time()
        self.electrodes = [[] for j in range(self.virtual_number_axons)]
        if self.number_contact_points == 1:
            for k in range(self.virtual_number_axons):
                filename = "electrode_"+str(k)+".dat"
                self.electrodes[k] = np.loadtxt(directory + filename, unpack=True)
                #for i in range(self.number_contact_points*len(self.recording_elec_pos)):
                    #for j in range(self.number_elecs):
                        #self.electrodes[k].append(np.loadtxt(directory + filename, unpack=True, usecols=[i*self.number_elecs+j]))
        else:
            for k in range(self.virtual_number_axons):
                filename = "electrode_"+str(k)+".dat"
                self.electrodes[k] = np.loadtxt(directory + filename, unpack=True)
                        
        print "loaded in: "+str(time.time()-t0)








    def compute_CAP1D_fromfile(self):
        self.load_electrodes()
        self.CAP = []
        self.sum_CAP = 0
        for i in range(0,len(self.recording_elec_pos)*self.number_contact_points,len(self.recording_elec_pos)):
            if len(self.recording_elec_pos) == 1:
                self.CAP.append(self.electrodes[0][i])
            elif len(self.recording_elec_pos) == 2:
                self.CAP.append(self.electrodes[0][i]-self.electrodes[0][i+1])
        for i in range(len(self.recording_elec_pos),len(self.recording_elec_pos)*self.number_contact_points,len(self.recording_elec_pos)):
            for j in range(1,self.virtual_number_axons):
                if len(self.recording_elec_pos) == 1:
                    self.CAP[int(float(i))] +=  self.electrodes[j][i]
                elif len(self.recording_elec_pos) == 2:
                    self.CAP[int(float(i)/2)] += (self.electrodes[j][i]- self.electrodes[j][i+1])
                    
            if len(self.recording_elec_pos) == 1:
                self.sum_CAP += self.CAP[int(float(i))]
            elif len(self.recording_elec_pos) == 2:                
                self.sum_CAP += self.CAP[int(float(i)/2)]

                        

    def compute_CAP2D_fromfile(self):
        directory = getDirectoryName("elec", **self.saveParams)

        temp = time.time()
        CAP = []
        print "loading electrode"
        t0 = time.time()

        filename = "electrode_"+str(0)+".dat"
        electrodesData = np.loadtxt(directory + filename, unpack=True)
        print "loaded in: "+str(time.time()-t0)
        self.sum_CAP = np.zeros((self.number_elecs,len(electrodesData[0])))
        for i in range(self.number_contact_points):
            for j in range(self.number_elecs):
                CAP.append(electrodesData[i*self.number_elecs+j])

        del electrodesData
        for k in range(1,self.virtual_number_axons):
            t0 = time.time()
            filename = "electrode_"+str(k)+".dat"
            electrodesData = np.loadtxt(directory + filename, unpack=True)
            print "electrode_"+str(k)+ "loaded in: "+str(time.time()-t0)
            for j in range(self.number_elecs):
                for i in range(self.number_contact_points):
                    CAP[i*self.number_elecs+j] += electrodesData[i*self.number_elecs+j]
                    self.sum_CAP[j,:] += CAP[i*self.number_elecs+j]
            del electrodesData

        elapsed = time.time()-temp
        print "Elapsed time to compute CAP:" + str(elapsed)


    def setup_recording_elec(self):
        if (self.number_elecs == 1):
            if len(self.recording_elec_pos) == 1:
                X = np.zeros(self.number_contact_points)+self.recording_elec_pos
            elif len(self.recording_elec_pos) == 2:
                X = np.repeat(self.recording_elec_pos,self.number_contact_points,axis=0)
            elif len(self.recording_elec_pos) > 2 or len(self.recording_elec_pos) == 0:
                raise TypeError("Only monopolar and bipolar recording are supported")                
        else:
            if len(self.recording_elec_pos) >1:
                raise TypeError("Please use only the monopolar configuration of 'recording_elec' to record along the bundle") 
            X1 = [np.linspace(0, self.recording_elec_pos[0], self.number_elecs)] 
            X = np.repeat(X1,self.number_contact_points, axis=0)
            X = X.flatten()
        angles = np.linspace(0,360, self.number_contact_points, endpoint = False)
        Y1 = np.round(self.radius_bundle*np.cos(angles*np.pi/180.0),2)
        Z1 = np.round(self.radius_bundle*np.sin(angles*np.pi/180.0),2)
        Y = np.repeat(Y1,self.number_elecs*len(self.recording_elec_pos))
        Z = np.repeat(Z1,self.number_elecs*len(self.recording_elec_pos))
        N = np.empty((self.number_contact_points*self.number_elecs*len(self.recording_elec_pos), 3))
        for i in xrange(N.shape[0]):
            N[i,] = [1, 0, 0] #normal vec. of contacts
        return [angles,X,Y,Z,N]
        
    def build_disk(self,number_of_axons,radius_bundle):
        ### AXONS POSITIONS ON THE DISK ###
        """
        Partly from http://blog.marmakoide.org/?p=1
        """
        n = number_of_axons 

        radius = radius_bundle* np.sqrt(np.arange(n) / float(n))

        golden_angle = np.pi * (3 - np.sqrt(5))
        theta = golden_angle * np.arange(n)

        self.axons_pos = np.zeros((n, 2))
        self.axons_pos[:,0] = np.cos(theta)
        self.axons_pos[:,1] = np.sin(theta)
        self.axons_pos *= radius.reshape((n, 1))
        
    def get_filename(self, recordingType):

        directory = getDirectoryName(recordingType, **self.saveParams)
        if not os.path.exists(directory):
            os.makedirs(directory)

        filename = 'recording.dat'

        number = 0
        while os.path.isfile(directory+filename):
            number += 1
            print "Be careful this file name already exist ! We concatenate a number to the name to avoid erasing your previous file."
            filename = str(number) + filename

        self.filename = directory+filename

        return self.filename

    def get_filename_for_draw(self):
        self.filename = "p_A"+str(self.p_A)+"_p_C"+str(self.p_C)+'nb_axons'+str(self.number_of_axons)+'bundle_radius'+str(self.radius_bundle)+'.dat'

        return self.filename

    def getDiam(self, axonType):

        if axonType == 'm':
            givenDiameter = self.myelinated_A['fiberD']

            if isinstance(givenDiameter, float) or isinstance(givenDiameter, int):
                return givenDiameter
            elif isinstance(givenDiameter, dict):
                # get diameter distribution
                fiberD = self.myelinated_A['fiberD']
                densities = fiberD['densities']

                # normalize it
                sum_d = float(sum(densities))
                normalize_densities = [x / sum_d for x in densities]

                # draw one diameter value from the distribution
                draw_diam = np.random.choice(len(normalize_densities),1,p = normalize_densities)

                # why add 2.8?
                axonD = fiberD['diameters'][draw_diam]+2.8

                # choose the closest from existing axonD
                axonD_choices = [3.4,4.6,6.9,8.1,9.2,10.4,11.5,12.7] # assuming the user consider usually in the axon diameter
                diff_axonD = [abs(x-axonD) for x in axonD_choices]
                fiberD_choices = [5.7, 7.3, 8.7, 10.0, 11.5, 12.8, 14.0, 15.0, 16.0]
                fiberD = fiberD_choices[np.argmin(diff_axonD)]
                diam = fiberD
            else:
                raise('Something is wrong with your given axon diameter for myelinated axons.')

        elif  axonType == 'u':
            givenDiameter = self.unmyelinated['diam']

            if isinstance(givenDiameter, float) or isinstance(givenDiameter, int):
                return givenDiameter
            elif isinstance(givenDiameter, dict):
                fiberD = self.unmyelinated['diam']
                densities = fiberD['densities']

                sum_d = float(sum(densities))
                normalize_densities = [x / sum_d for x in densities]

                draw_diam = np.random.choice(len(normalize_densities),1,p = normalize_densities)
                D = fiberD['diameters'][draw_diam]
                diam = D
            else:
                raise('Something is wrong with your given axon diameter for unmyelinated axons.')


        else:
            raise('Wrong axon type given to function getDiam. Valid ones: u or m')

        return diam

    # def draw_distrib(self):
    #     startTime = time.time()
    #
    #     draw = np.random.choice(2,self.number_of_axons,p = [self.p_A, self.p_C])
    #     diams = []
    #     for i in range(len(draw)):
    #         if (isinstance(self.myelinated_A['fiberD'],float) and (isinstance(self.unmyelinated['diam'],float) or isinstance(self.unmyelinated['diam'],int))):
    #             pass
    #         else:
    #             if draw[i] == 0:
    #                 fiberD = self.myelinated_A['fiberD'];
    #                 densities = fiberD['densities'];
    #                 sum_d = float(sum(densities))
    #                 normalize_densities = [x / sum_d for x in densities]
    #                 draw_diam = np.random.choice(len(normalize_densities),1,p = normalize_densities)
    #                 axonD = fiberD['diameters'][draw_diam]+2.8
    #                 # choose the closest from existing axonD
    #                 axonD_choices = [3.4,4.6,6.9,8.1,9.2,10.4,11.5,12.7] # assuming the user consider usually in the axon diameter
    #                 diff_axonD = [abs(x-axonD) for x in axonD_choices]
    #                 fiberD_choices = [5.7, 7.3, 8.7, 10.0, 11.5, 12.8, 14.0, 15.0, 16.0]
    #                 fiberD = fiberD_choices[np.argmin(diff_axonD)]
    #                 diams.append(fiberD)
    #             else:
    #                 """print "Unmyelinated in the function"
    #                 print self.unmyelinated
    #                 print "\n" """
    #                 diam = self.unmyelinated['diam']
    #                 densities = diam['densities'];
    #                 sum_d = float(sum(densities))
    #                 normalize_densities = [x / sum_d for x in densities]
    #                 draw_diam = np.random.choice(len(normalize_densities),1,p = normalize_densities)
    #                 D = diam['diameters'][draw_diam]
    #                 diams.append(D)
    #
    #     print 'Time to draw from diameter distribution: ' + str(time.time() - startTime)
    #
    #
    #     parameters_to_save = {'p_A': self.p_A, 'p_C': self.p_C, 'radius_bundle': self.radius_bundle, 'number_of_axons':self.number_of_axons, 'unmyelinated': self.unmyelinated, 'myelinated': self.myelinated_A}
    #     header = repr(parameters_to_save)
    #     filename = self.get_filename_for_draw()
    #     directory = "FOR_PAPER/CAP2D/draws/"
    #     if not os.path.exists(directory):
    #         os.makedirs(directory)
    #     DataOut = np.column_stack( (draw, diams))
    #     filename_org = filename
    #     number  = 0
    #     while os.path.isfile(directory+filename):
    #         number += 1
    #         print "Be careful this file name already exist ! We concatenate a number to the name to avoid erasing your previous file."
    #         filename = str(number) + filename_org
    #     np.savetxt(directory+filename, DataOut, header=header)
    #     return [draw, diams]


    def load_distrib(self):
        #directory = "draws/biphasic/"
        # saveParams={'elecCount': len(self.recording_elec_pos), 'dt': h.dt, 'tStop': h.tstop, 'p_A': self.p_A,
        #             'myelinatedDiam': self.myelinated_A['fiberD'], 'unmyelinatedDiam': self.unmyelinated['diam'],
        #             'L': self.unmyelinated['L'], 'stimType': self.stim_type, 'stimWaveform' : self.waveform,
        #             'stimDutyCycle': self.duty_cycle, 'stimAmplitude' : self.amp}
        directory = getDirectoryName("draw", **self.saveParams)
        filename = self.get_filename_for_draw()
        draw = np.loadtxt(directory +filename, unpack=True, usecols=[0])
        diams = np.loadtxt(directory +filename, unpack=True, usecols=[1])
        return [draw,diams]


## METHOD USED TO SET SOME HEADER PARAMETERS ###
def fiberD_dependent_param(fiberD, nodelength, paralength1):
    if (fiberD==5.7):            
        deltax=500
        paralength2=35            
    if (fiberD==7.3):
        deltax=750
        paralength2=38
    if (fiberD==8.7):
        deltax=1000
        paralength2=40
    if (fiberD==10.0):
        deltax=1150
        paralength2=46
    if (fiberD==11.5):
        deltax=1250
        paralength2=50
    if (fiberD==12.8):
        deltax=1350
        paralength2=54
    if (fiberD==14.0):
        deltax=1400
        paralength2=56
    if (fiberD==15.0):
        deltax=1450
        paralength2=58
    if (fiberD==16.0):
        deltax=1500
        paralength2=60
        
    interlength = (deltax-nodelength-(2*paralength1)-(2*paralength2))/6.0
    return [paralength2,interlength]
