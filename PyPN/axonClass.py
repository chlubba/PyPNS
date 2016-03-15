import neuron
from neuron import h
import LFPy
import numpy as np # for arrays managing
import math
import os

import createGeometry


# Some methods in the Axon class are based on existing methods in the Python package LFPY

class Axon(object):
    # Own initialization method of the superclass Axon
    # rec_v: set voltage recorders True or False
    # layout3D: either "DEFINE_SHAPE" or "PT3D" using hoc corresponding function


    def __init__(self, layout3D, rec_v, name, fiberD, coord, temperature, tStop, timeRes):
        self.layout3D = layout3D
        self.rec_v = rec_v
        self.name = name
        self.fiberD = fiberD
        self.coord = coord
        self.synapse = []
        # LFPy initilizations
        self.verbose = False
        self.dotprodresults = None # Not in class cell of LFPY but present in run_simulation as cell.dotprodresults
        self.exMechVars = []

        # params for NEURON simulation
        self.temperature = temperature # set temperature in celsius
        self.tStop = tStop # set simulation duration (ms)
        self.timeRes = timeRes # set time step (ms)

        self.create_cell_props_for_LFPy(tStop, timeRes)

    def calc_totnsegs(self):
        # Calculate the number of segments in the allseclist
        i = 0
        for sec in self.allseclist:
            i += sec.nseg

        return i

    # this function and the following one are only there for compliance with the LFPy package
    def create_cell_props_for_LFPy(self, tStop, timeRes):
        self.tstartms = 0
        self.tstopms = tStop
        self.timeres_NEURON = timeRes
        self.timeres_python = timeRes
    def _loadspikes(self):
        pass

    def append_ex_mech_vars(self, exMechVars):
        self.exMechVars.append(exMechVars)

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
        # Loop over allseclist to determine area, diam, xyz-start- and
        # endpoints, embed geometry to self object
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

        # Record membrane currents for all segments

        self.memireclist = neuron.h.List()
        for sec in self.allseclist:
            for seg in sec:
                memirec = h.Vector(int(h.tstop/h.dt +1))
                memirec.record(seg._ref_i_membrane)
                self.memireclist.append(memirec)
    def set_voltage_recorders(self):

        # Record voltage for all segments (not from LFPy sources)

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

        # Fetch the vectors from the memireclist and calculate self.imem
        # containing all the membrane currents.

        self.imem = np.array(self.memireclist)
        for i in range(self.imem.shape[0]):
            self.imem[i, ] *= self.area[i] * 1E-2
        self.memireclist = None
        del self.memireclist

    def collect_tvec(self):

        # Set the tvec to be a monotonically increasing numpy array after sim.
        self.tvec = np.arange(h.tstop /h.dt + 1)*h.dt

    def simulate(self, electrode=None, rec_imem=True,
                 rec_variables=[], variable_dt=False, atol=0.001,
                 to_memory=True, to_file=False, file_name=None,
                 dotprodcoeffs=None):


        # This is the main function running the simulation of the NEURON model.
        # Start NEURON simulation and record variables specified by arguments.
        #
        # Arguments:
        # ::
        #
        #     electrode:  Either an LFPy.RecExtElectrode object or a list of such.
        #                 If supplied, LFPs will be calculated at every time step
        #                 and accessible as electrode.LFP. If a list of objects
        #                 is given, accessible as electrode[0].LFP etc.
        #     rec_imem:   If true, segment membrane currents will be recorded
        #                 If no electrode argument is given, it is necessary to
        #                 set rec_imem=True in order to calculate LFP later on.
        #                 Units of (nA).
        #     rec_v:      record segment voltages (mV)
        #     rec_istim:  record currents of StimIntraElectrode (nA)
        #     rec_variables: list of variables to record, i.e arg=['cai', ]
        #     variable_dt: boolean, using variable timestep in NEURON
        #     atol:       absolute tolerance used with NEURON variable timestep
        #     to_memory:  only valid with electrode, store lfp in -> electrode.LFP
        #     to_file:    only valid with electrode, save LFPs in hdf5 file format
        #     file_name:  name of hdf5 file, '.h5' is appended if it doesnt exist
        #     dotprodcoeffs :  list of N x Nseg np.ndarray. These arrays will at
        #                 every timestep be multiplied by the membrane currents.
        #                 Presumably useful for memory efficient csd or lfp calcs



        if rec_imem:
            self.set_imem_recorders()
        if self.rec_v:
            self.set_voltage_recorders()

        h.celsius = self.temperature # set temperature in celsius
        LFPy.run_simulation._run_simulation_with_electrode(self, electrode, variable_dt, atol,
                                                   to_memory, to_file, file_name,
                                                   dotprodcoeffs)
        if rec_imem:
            self.calc_imem()


    # Equivalent methods interpxyz and setrx from the xtra mechanism available on the NEURON website from Ted Carnevale
    # Setrx has been modified to integrate the use of multipolar electrodes
    def interpxyz(self):
        # interpolated data, spaced at regular intervals

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


class Unmyelinated(Axon):

    # name: axon name (for neuron)
    #
    # nsegs_method: ['lambda100']/'lambda_f'/'fixed_length': nseg rule
    # max_nsegs_length: [None]: max segment length for method 'fixed_length'
    # lambda_f: [100]: AC frequency for method 'lambda_f'
    # d_lambda: [0.1]: parameter for d_lambda rule
    #
    # diam: Axon diameter (micrometer)
    # cm : Specific membrane capacitance (microfarad/cm2)
    # Ra: Specific axial resistance (Ohm cm)
    # coord: y,z coordinates to spatially place the axon once created
    # layout3D: either "DEFINE_SHAPE" or "PT3D" using hoc corresponding function
    # rec_v: set voltage recorders True or False

    def __init__(self, fiberD, coord, tStop, timeRes, temperature=33, cm=1.0, Ra=200.0, name="unmyelinated_axon", layout3D="PT3D", rec_v=True, hhDraw=False, nsegs_method='lambda100', lambda_f=100, d_lambda=0.1, max_nsegs_length=None):
        super(Unmyelinated,self).__init__(layout3D, rec_v, name, fiberD, coord, temperature, tStop, timeRes)

        self.L = createGeometry.length_from_coords(coord)
        self.cm = cm
        self.Ra = Ra
        self.hhDraw = hhDraw
        self.nsegs_method = nsegs_method
        self.lambda_f = lambda_f
        self.d_lambda = d_lambda
        self.max_nsegs_length = max_nsegs_length

        self.v_init = -64.975

        print "Unmyelinated axon diameter: " + str(self.fiberD)

    def create_neuron_object(self):
        self.axon = h.Section(name = str(self.name))
        self.allseclist = h.SectionList()
        self.allseclist.append(sec = self.axon)

        self.axon.insert('extracellular')
        self.axon.insert('xtra')
        self.axon_update_property()

        self.set_nsegs(self.nsegs_method, self.lambda_f, self.d_lambda, self.max_nsegs_length)

        if (self.layout3D == "DEFINE_SHAPE"):
            h.define_shape()
        elif (self.layout3D == "PT3D"):

            if True:
                h.pt3dclear(sec=self.axon)
                for i in range(np.shape(self.coord)[0]):
                    h.pt3dadd(self.coord[i,0], self.coord[i,1], self.coord[i,2], self.fiberD, sec=self.axon)
            else:
                h.pt3dadd(0, self.coord[0], self.coord[1], self.fiberD, sec=self.axon)
                h.pt3dadd(self.L, self.coord[0], self.coord[1], self.fiberD, sec=self.axon)
        else:
            raise NameError('layout3D only "DEFINE_SHAPE" or "PT3D"')

        for sec_id in self.allseclist:
            for seg in sec_id:
                h.setpointer(seg._ref_i_membrane, 'im', seg.xtra)
                h.setpointer(seg._ref_e_extracellular, 'ex', seg.xtra)
        self.interpxyz()
        self.collect_geometry()
        self.calc_midpoints()
        self.channel_init() # default values of hh channel are used

    def delete_neuron_object(self):

        for sec in self.allseclist:
            for seg in sec:
                seg = None
            sec = None
        self.allseclist = None
        self.axon = None

        if not self.synapse == []:
            self.synapse = None
            self.vecStim = None
            self.netCon = None

        try:
            for memirec in self.memireclist:
                memirec = None
            self.memireclist = None
        except:
            pass

        try:
            for vrec in self.vreclist:
                vrec = None
            self.vreclist = None
        except:
            pass

        # also delete unnecessary data that will no longer be used to keep the pickled file small
        self.imem = None

        # neuron objects inserted by excitation mechanisms
        for exMechVars in self.exMechVars:
            for i in range(len(exMechVars)):
                exMechVars[i] = None



    def axon_update_property(self):
        self.axon.L = self.L
        self.axon.diam = self.fiberD
        self.axon.cm = self.cm
        self.axon.Ra = self.Ra





    def channel_init(self, gna=0.120, gk=0.036, gl=0.0003, ena=50, ek=-77, el=-54.3):

        # g unit S/cm2 and e unit mV
        # gna: maximum specific sodium channel conductance
        # gk: maximum specific potassium channel conductance
        # gl: maximum specific leakage conductance
        # ena: reversal potential for the sodium channel
        # ek: reversal potential for the potassium channel
        # el: reversal potential for the leakage channel

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
        # Set number of segments per section according to the lambda-rule,
        # or according to maximum length of segments
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
        # Set the number of segments for section according to the
        # d_lambda-rule for a given input frequency
        #     frequency: float, frequency at whihc AC length constant is computed
        #     d_lambda: float,

        for sec in self.allseclist:
            sec.nseg = int((sec.L / (d_lambda*h.lambda_f(frequency,
                                                           sec=sec)) + .9)
                / 2 )*2 + 1
            print "Number of segments for unmyelinated axon via d_lambda: "+ str(sec.nseg)
        if self.verbose:
            print(("set nsegs using lambda-rule with frequency %i." % frequency))

    def set_nsegs_lambda100(self, d_lambda=0.1):
        # Set the numbers of segments using d_lambda(100)
        self.set_nsegs_lambda_f(frequency=100, d_lambda=d_lambda)

    def set_nsegs_fixed_length(self, maxlength):
        # Set nseg for sections so that every segment L < maxlength
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
    def __init__(self, fiberD, coord, tStop, timeRes, temperature=40, name="myelinated_axonA", layout3D="PT3D", rec_v=True):
        super(Myelinated,self).__init__(layout3D, rec_v, name, fiberD, coord, temperature, tStop, timeRes)

        self.v_init = -80

        print 'Myelinated fiber diameter: ' + str(self.fiberD)

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


        # length from the middle of one node to the middle of the next
        self.lengthOneCycle = self.nodelength + self.interlength*6 + 2*self.paralength1 + 2*self.paralength2

        # length of the whole axon
        self.L = createGeometry.length_from_coords(coord)

        # number of nodes
        self.axonnodes = int(math.ceil(self.L/self.lengthOneCycle))

        # number of remaining sections
        self.paranodes1= 2*(self.axonnodes-1)
        self.paranodes2= 2*(self.axonnodes-1)
        self.axoninter= 6*(self.axonnodes-1)
        self.axontotal= self.axonnodes+self.paranodes1+self.paranodes2+self.axoninter

    def createSingleNode(self, nodeType):

        if nodeType == 'n':
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

        if nodeType == 'm':
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

        if nodeType == 'f':
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

        if nodeType == 's':
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

    def create_neuron_object(self):

        self.allseclist = h.SectionList()

        self.nodes = []
        self.MYSAs = []
        self.FLUTs = []
        self.STINs = []

        # iterate through nodes in the way they are ordered in the axon in order to have an ordered allsectionlist
        nodeSequence = ['n', 'm', 'f', 's', 's', 's', 's', 's', 's', 'f', 'm']
        for i in range(self.axonnodes-1):
            for j in range(len(nodeSequence)):
                nodeType = nodeSequence[j]
                self.createSingleNode(nodeType)
        self.createSingleNode('n')

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

            lengthArray = np.concatenate(([self.nodelength, self.paralength1, self.paralength2], np.multiply(np.ones(6), self.interlength), [self.paralength2, self.paralength1]))

            # get sections in order (allseclist not subscriptable since NEURON object. But then, iterable.)
            sectionArray = []
            for sec in self.allseclist:
                sectionArray.append(sec)


            # DO-WHILE emulation
            sectionIndex = 0
            coordCounter = 0

            # get the segment from the axon coordinates defined by createRandomAxon
            direction = self.coord[coordCounter + 1, :] - self.coord[coordCounter, :]

            # calculate the length and the normed direction vector
            lengthGuideSegment = np.linalg.norm(direction)
            directionNorm = direction / np.linalg.norm(direction)

            # while on the same guide segment
            sectionTypeIndex = sectionIndex % 11
            sectionLength = lengthArray[sectionTypeIndex]

            # variables to keep track of remaining section lengths and guide segment lengths
            cumulatedLengthOnGuideSeg = 0
            cumulatedLengthFromSection = 0

            # set startpoint of axon /first axon section
            section = sectionArray[sectionIndex]
            h.pt3dclear(sec=section)
            coord=self.coord[0,:]
            h.pt3dadd(coord[0], coord[1], coord[2], section.diam, sec=section)

            # print 'Section ' + str(sectionIndex) + ' started at coords ' + str(coord) + '.'
            lengthReached = 0

            while sectionIndex < self.axontotal:

                # if the axon guide segment is longer than the remaining axon section, go to next section
                if cumulatedLengthOnGuideSeg + (sectionLength - cumulatedLengthFromSection) < lengthGuideSegment:

                    # set coordinates for section as lying in axon guide
                    # coord = self.coord[coordCounter,:] + directionNorm * (sectionLength - cumulatedLengthFromSection)
                    coord = coord + directionNorm * (sectionLength - cumulatedLengthFromSection)
                    lengthReached += (sectionLength - cumulatedLengthFromSection)

                    # set endpoint of section
                    h.pt3dadd(coord[0], coord[1], coord[2], section.diam, sec=section)
                    # print 'Section ' + str(sectionIndex) + ' ended at coords' + str(coord) + '.'
                    # print 'Direction normed : ' + str(directionNorm) + '\n'

                    # update how far we got on the segment of the guide
                    cumulatedLengthOnGuideSeg += (sectionLength - cumulatedLengthFromSection)

                    # step to next segment
                    sectionIndex += 1

                    if sectionIndex >= self.axontotal:
                        break

                    # get the length
                    sectionTypeIndex = sectionIndex % 11
                    sectionLength = lengthArray[sectionTypeIndex]

                    # set startpoint of next section
                    section = sectionArray[sectionIndex]
                    h.pt3dclear(sec=section)
                    h.pt3dadd(coord[0], coord[1], coord[2], section.diam, sec=section)
                    # print 'Section ' + str(sectionIndex) + ' started at coords' + str(coord) + '.'
                    # print 'Length reached : '+str(lengthReached)

                    section.connect(sectionArray[sectionIndex-1],1,0)

                    cumulatedLengthFromSection  = 0

                # if section does not fit into remaining axon guide segment, go to next axon guide segment
                else:

                    # substract the remaining length on the former segment from the section length to see what is left
                    # for consecutive segments
                    cumulatedLengthFromSection += lengthGuideSegment - cumulatedLengthOnGuideSeg

                    lengthReached += np.linalg.norm(self.coord[coordCounter + 1,:] - coord)

                    # get coords of node between guide segments
                    coord = self.coord[coordCounter + 1,:]

                    h.pt3dadd(coord[0], coord[1], coord[2], section.diam, sec=section)
                    # print 'Section ' + str(sectionIndex) + ' has additional corner at coords' + str(coord) + '.'
                    # print 'Length reached : '+str(lengthReached)

                    # check out next guide segment
                    coordCounter += 1

                    # calculate the next normed direction vector
                    direction = self.coord[coordCounter + 1, :] - self.coord[coordCounter, :]
                    directionNorm = direction / np.linalg.norm(direction)
                    lengthGuideSegment = np.linalg.norm(direction)

                    # new guide segment -> no length consumed yet.
                    cumulatedLengthOnGuideSeg = 0

        else:
            raise NameError('layout3D only "DEFINE_SHAPE" or "PT3D"')


        self.totnsegs = self.calc_totnsegs()
        for sec_id in self.allseclist:
            for seg in sec_id:
                h.setpointer(seg._ref_i_membrane, 'im', seg.xtra)
                h.setpointer(seg._ref_e_extracellular, 'ex', seg.xtra)
        self.interpxyz()
        self.collect_geometry()
        self.calc_midpoints()

    def delete_neuron_object(self):
        for sec in self.allseclist:
            for seg in sec:
                seg = None
            sec = None
        self.allseclist = None

        self.nodes = None
        self.FLUTs = None
        self.MYSAs = None
        self.STINs = None

        if not self.synapse == []:
            self.synapse = None
            self.vecStim = None
            self.netCon = None

        if not self.synapse == []:
            self.synapse = None
            self.vecStim = None
            self.netCon = None

        try:
            for memirec in self.memireclist:
                memirec = None
            self.memireclist = None
        except:
            pass

        try:
            for vrec in self.vreclist:
                vrec = None
            self.vreclist = None
        except:
            pass

        # also delete unnecessary data that will no longer be used to keep the pickled file small
        self.imem = None

        # do the following better, soon!! have a separate excitationMechanism list in each axon

        for exMechVars in self.exMechVars:
            for i in range(len(exMechVars)):
                exMechVars[i] = None