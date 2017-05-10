import neuron
from neuron import h
err = h('load_file("noload.hoc")')
# if not err == 1:
#     raise ImportError('The needed mod-files were not found. Run "nrnivmodl" in the mod-directory and copy the result to the project directory.')
import numpy as np # for arrays managing
import math
import os
from scipy import interpolate
import numbers

import createGeometry


# Some methods in the Axon class are based on existing methods in the Python package LFPY

class Axon(object):

    def __init__(self, layout3D, rec_v, name, fiberD, coord, temperature, v_init, tStop, timeRes, numberOfSavedSegments):

        # properties of the axon
        self.name = name
        self.fiberD = fiberD
        self.coord = coord
        self.v_init = v_init

        # configuration params of the simulation
        self.layout3D = layout3D
        self.rec_v = rec_v
        self.verbose = False

        # excitation mechanisms added to the axon
        # self.synapse = []
        self.exMechVars = []

        # params for NEURON simulation
        self.temperature = temperature # set temperature in celsius
        self.tStop = tStop # set simulation duration (ms)
        self.timeRes = timeRes # set time step (ms)
        self.numberOfSavedSegments = numberOfSavedSegments

        # initialize variables for geometry
        self.xstart = None
        self.ystart = None
        self.zstart = None
        self.xend = None
        self.yend = None
        self.zend = None
        self.area = None
        self.diam = None
        self.length = None

        # initialize save variables for NEURON objects (still in NEURON format)
        self.memireclist = None
        self.vreclist = None
        self.allseclist = h.SectionList()
        # self.axon = None

        # save variables for numpy format, imem being scaled to nA (from mA/cm^2)
        self.imem = None
        self.tvec = None
        self.totnsegs = None

        # save the type of segments in an array as to identify them afterwards
        self.segmentTypes = None

    def calc_totnsegs(self):
        # Calculate the number of segments in the allseclist (only possible if NEURON simulation has run)
        i = 0
        for sec in self.allseclist:
            i += sec.nseg

        return i

    def append_ex_mech_vars(self, exMechVars):
        self.exMechVars.append(exMechVars)

    def collect_geometry(self):
        """
        Collects x, y, z-coordinates from NEURON

        Returns: Nothing, write into Axon-class
        """

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
        """
        Loop over allseclist to determine area, diam, xyz-start- and
        endpoints, embed geometry to self object

        Returns:

        """

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
        self.xmid = .5*(self.xstart+self.xend).flatten()
        self.ymid = .5*(self.ystart+self.yend).flatten()
        self.zmid = .5*(self.zstart+self.zend).flatten()

    def set_imem_recorders(self):

        # Record membrane currents for all segments (scale with area and a factor to adjust the unit)
        self.memireclist = neuron.h.List()
        for sec in self.allseclist:
            for seg in sec:
                memirec = h.Vector(int(h.tstop/h.dt +1))
                memirec.record(seg._ref_i_membrane)
                self.memireclist.append(memirec)

    def set_voltage_recorders(self):

        # Record voltage for all segments
        self.vreclist = h.List()

        for sec in self.allseclist:
            # address the problem of the high number of segments necessary to compute the accurate AP propagation
            # in the unmyelinated axon case by limiting the number of monitored segments.
            if sec.nseg > self.numberOfSavedSegments:
                for i in range(1,self.numberOfSavedSegments+1):
                    vrec = h.Vector(int(h.tstop/h.dt+1))
                    vrec.record(sec(float(i)/self.numberOfSavedSegments)._ref_v)
                    self.vreclist.append(vrec)
            else:
                for seg in sec:
                    vrec = h.Vector(int(h.tstop/h.dt+1))
                    vrec.record(seg._ref_v)
                    self.vreclist.append(vrec)

    def calc_imem(self):
        """
        Fetch the vectors from the memireclist and calculate self.imem
        containing all the membrane currents. From LFPy

        Returns: -

        """

        self.imem = np.array(self.memireclist)

        for i in range(self.imem.shape[0]):
            self.imem[i, ] *= self.area[i] * 1E-2 # * 1E2 #
        self.memireclist = None
        del self.memireclist

    def collect_tvec(self):

        # Set the tvec to be a monotonically increasing numpy array after sim.
        self.tvec = np.arange(h.tstop /h.dt + 1)*h.dt

    def set_nsegs(self):
        """
        Set number of segments per section according to the lambda-rule, or according to maximum length of segments

        """

        d_lambda = 0.1
        frequency = 100

        for sec in self.allseclist:
            sec.nseg = int((sec.L / (d_lambda*h.lambda_f(frequency, sec=sec)) + .9) / 2 ) * 2 + 1

        if self.verbose:
            print("set nsegs using lambda-rule with frequency %i." % frequency)

        self.totnsegs = self.calc_totnsegs()


    def position_sections_in_neuron(self):

        # get sections in order (allseclist not subscriptable since NEURON object. But then, iterable.)
        sectionArray = []
        for sec in self.allseclist:
            sectionArray.append(sec)
        numSecs = np.shape(sectionArray)[0]

        # DO-WHILE emulation
        sectionIndex = 0
        coordCounter = 0

        # get the segment from the axon coordinates defined by createRandomAxon
        direction = self.coord[coordCounter + 1, :] - self.coord[coordCounter, :]

        # calculate the length and the normed direction vector
        lengthGuideSegment = np.linalg.norm(direction)
        directionNorm = direction / np.linalg.norm(direction)

        # variables to keep track of remaining section lengths and guide segment lengths
        cumulatedLengthOnGuideSeg = 0
        cumulatedLengthFromSection = 0

        # set startpoint of axon /first axon section
        section = sectionArray[sectionIndex]

        # get length
        sectionLength = section.L

        # set actual coordinate values
        h.pt3dclear(sec=section)
        coord = self.coord[0, :]
        h.pt3dadd(coord[0], coord[1], coord[2], section.diam, sec=section)

        # print 'Section ' + str(sectionIndex) + ' started at coords ' + str(coord) + '.'
        lengthReached = 0

        while sectionIndex < numSecs:

            # if the axon guide segment is longer than the remaining axon section, go to next section
            if cumulatedLengthOnGuideSeg + (sectionLength - cumulatedLengthFromSection) < lengthGuideSegment:

                # set coordinates for section as lying in axon guide
                # coord = self.coord[coordCounter,:] + directionNorm * (sectionLength - cumulatedLengthFromSection)
                coord += directionNorm * (sectionLength - cumulatedLengthFromSection)
                lengthReached += (sectionLength - cumulatedLengthFromSection)

                # set endpoint of section
                h.pt3dadd(coord[0], coord[1], coord[2], section.diam, sec=section)

                # update how far we got on the segment of the guide
                cumulatedLengthOnGuideSeg += (sectionLength - cumulatedLengthFromSection)

                # step to next segment
                sectionIndex += 1

                if sectionIndex >= numSecs:
                    break

                # set startpoint of next section
                section = sectionArray[sectionIndex]

                # get length
                sectionLength = section.L

                # set actual coordinate values
                h.pt3dclear(sec=section)
                h.pt3dadd(coord[0], coord[1], coord[2], section.diam, sec=section)

                section.connect(sectionArray[sectionIndex - 1], 1, 0)

                cumulatedLengthFromSection = 0

            # if section does not fit into remaining axon guide segment, go to next axon guide segment
            else:

                # substract the remaining length on the former segment from the section length to see what is left
                # for consecutive segments
                cumulatedLengthFromSection += lengthGuideSegment - cumulatedLengthOnGuideSeg

                lengthReached += np.linalg.norm(self.coord[coordCounter + 1, :] - coord)

                # get coords of node between guide segments
                coord = self.coord[coordCounter + 1, :]

                h.pt3dadd(coord[0], coord[1], coord[2], section.diam, sec=section)

                # check out next guide segment
                if coordCounter < np.shape(self.coord)[0]-2:
                    coordCounter += 1
                # else:
                #     print 'whoops, axon guide not long enough for axon.'

                # calculate the next normed direction vector
                direction = self.coord[coordCounter + 1, :] - self.coord[coordCounter, :]
                directionNorm = direction / np.linalg.norm(direction)
                lengthGuideSegment = np.linalg.norm(direction)

                # new guide segment -> no length consumed yet.
                cumulatedLengthOnGuideSeg = 0

    def create_neuron_object(self):


        # # general things to do
        # self.allseclist = h.SectionList()
        # self.allseclist.append(sec = self.axon)

        # calculate number of segments per section
        self.set_nsegs()

        for sec in self.allseclist:

            # insert extracellular mechanism for imem
            sec.insert('extracellular')
            # insert xtra mechanism for extracellular stimulation
            sec.insert('xtra')

            for seg in sec:
                h.setpointer(seg._ref_i_membrane, 'im', seg.xtra)
                h.setpointer(seg._ref_e_extracellular, 'ex', seg.xtra)

        # get geometry from NEURON object
        self.interpxyz()
        self.collect_geometry()
        self.calc_midpoints()

    def delete_neuron_object(self):

        for sec in self.allseclist:
            for seg in sec:
                seg = None
            sec = None
        self.allseclist = None

        # if not self.synapse == []:
        #     self.synapse = None
        #     self.vecStim = None
        #     self.netCon = None

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

    def simulate(self, rec_imem=True):

        # before LFPy was used to simulate, this is not necessary.

        # self, electrode = None, rec_imem = True,
        # rec_variables = [], variable_dt = False, atol = 0.001,
        # to_memory = True, to_file = False, file_name = None,
        # dotprodcoeffs = None):

        # This is the main function running the simulation of the NEURON model.
        # Start NEURON simulation and record variables specified by arguments.
        #
        # Arguments:
        # ::
        #
        #     rec_imem:   If true, segment membrane currents will be recorded
        #                 If no electrode argument is given, it is necessary to
        #                 set rec_imem=True in order to calculate LFP later on.
        #                 Units of (nA).



        if rec_imem:
            self.set_imem_recorders()
        if self.rec_v:
            self.set_voltage_recorders()

        h.celsius = self.temperature # set temperature in celsius

        # LFPy.run_simulation._run_simulation_with_electrode(self, electrode, variable_dt, atol, to_memory, to_file, file_name, dotprodcoeffs)
        h.v_init = self.v_init
        h.finitialize(self.v_init)
        h.tstop = self.tStop

        # variable time step is indicated by a string as the timeRes parameter of the bundle/axon. Normally 'variable'
        if isinstance(self.timeRes, numbers.Number):
            h.dt = self.timeRes
        else:
            h('cvode_active(1)')
            h('cvode.atol(0.001)')

        h.run()

        if rec_imem:
            self.calc_imem()

    def interpxyz(self):
        """
        Equivalent methods interpxyz and setrx from the xtra mechanism available on the NEURON website from Ted Carnevale
        Setrx has been modified to integrate the use of multipolar electrodes

        interpolated data, spaced at regular intervals

        Returns:

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

    @staticmethod
    def setrx(stim_elec, rho=500, bipolar = False):
        """
        This function assumes a homogeneous outer medium for stimulation and sets the transmission resistance rx
        accordingly. While for the CNS this might be ok, it certainly isn't here in the PNS, so this function should
        not be used and could as well be removed.

        Args:
            stim_elec:
            axon_pos:
            rho:
            bipolar:

        Returns:

        """

        numPoints = stim_elec.shape[0]

        r = np.zeros(len(stim_elec))

        segCounter = 0
        # now expects xyc coords as arguments
        for sec in h.allsec():

            if h.ismembrane('xtra', sec=sec):

                for segInd, seg in enumerate(sec):

                    for j in range(numPoints):

                        [xe, ye, ze] = stim_elec[j,:]

                        # avoid nodes at 0 and 1 ends, so as not to override values at internal nodes
                        r[j] = math.sqrt(
                            math.pow(seg.x_xtra - xe, 2) + math.pow(seg.y_xtra - ye, 2) + math.pow(
                                seg.z_xtra - ze, 2))

                        # 0.01 converts rho's cm to um and ohm to megohm
                        # if electrode is exactly at a node, r will be 0
                        # this would be meaningless since the location would be inside the cell
                        # so force r to be at least as big as local radius

                        if r[j] == 0:
                            r[j] = seg.diam / 2.0

                    sum_r = 0
                    pointIndex = 0
                    while pointIndex < numPoints:
                        if bipolar:
                            sign = (-1)**pointIndex
                        else:
                            sign = 1
                        sum_r += sign / r[pointIndex]
                        pointIndex += 1

                    # seg.xtra.rx = (sec.Ra / 4.0 / math.pi) * sum_r * 0.01
                    seg.xtra.rx = (rho / 4.0 / math.pi) * sum_r * 0.01

                    # print "seg no %i, rx %10.8f" % (segCounter, (rho / 4.0 / math.pi) * sum_r * 0.01)
                    segCounter += 1

class Unmyelinated(Axon):

    """
    name: axon name (for neuron)
    nsegs_method: ['lambda100']/'lambda_f'/'fixed_length': nseg rule
    max_nsegs_length: [None]: max segment length for method 'fixed_length'
    lambda_f: [100]: AC frequency for method 'lambda_f'
    d_lambda: [0.1]: parameter for d_lambda rule

    diam: Axon diameter (micrometer)
    cm : Specific membrane capacitance (microfarad/cm2)
    Ra: Specific axial resistance (Ohm cm)
    coord: y,z coordinates to spatially place the axon once created
    layout3D: either "DEFINE_SHAPE" or "PT3D" using hoc corresponding function
    rec_v: set voltage recorders True or False
    """

    def __init__(self, fiberD, coord, tStop, timeRes, numberOfSavedSegments, temperature=33, v_init=-64.975, cm=1.0, Ra=200.0, name="unmyelinated_axon", layout3D="PT3D", rec_v=True, hhDraw=False): # , nsegs_method='lambda100', lambda_f=100, d_lambda=0.1, max_nsegs_length=None):
        super(Unmyelinated,self).__init__(layout3D, rec_v, name, fiberD, coord, temperature, v_init, tStop, timeRes, numberOfSavedSegments)

        self.L = createGeometry.length_from_coords(coord)
        self.cm = cm
        self.Ra = Ra
        self.hhDraw = hhDraw

        print "Unmyelinated axon diameter: " + str(self.fiberD)

        print 'Number of segments for unmyelinated axon: %i' % self.get_number_of_segs()

    def get_number_of_segs(self, d_lambda=0.1, lambda_freq=100):

        def lambda_f(freq, diam, Ra, cm):
            return 1e5 * np.sqrt(diam / (4 * np.pi * freq * Ra * cm))

        def approx_nseg_d_lambda(axon):
            return int((axon.L / (d_lambda * lambda_f(lambda_freq, axon.fiberD, axon.Ra, axon.cm)) + 0.9) / 2) * 2 + 1

        return approx_nseg_d_lambda(self)

    def set_segment_types(self):

        self.segmentTypes = np.tile('n', self.totnsegs)

    def create_neuron_object(self):

        # define the section properties of the unmyelinated axon
        self.axon = h.Section(name = str(self.name))
        self.axon.L = self.L
        self.axon.diam = self.fiberD
        self.axon.cm = self.cm
        self.axon.Ra = self.Ra

        # add section(s) to the section list
        self.allseclist.append(sec = self.axon)

        # transfer the geometry into the NEURON axon model
        self.position_sections_in_neuron()

        # do everything (number of segs per sec, channels, geometry from NEURON)
        super(Unmyelinated, self).create_neuron_object()

        # as number of segments per section are set dynamically depending on diameter, set the association here for
        # further processing
        self.set_segment_types()

        # add channels specific to unmyelinated axon
        self.channel_init() # default values of hh channel are used

    def delete_neuron_object(self):

        super(Unmyelinated, self).delete_neuron_object()
        self.axon = None

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





'''/*--------------------------------------------------------------------
2/02
Cameron C. McIntyre
SIMULATION OF PNS MYELINATED AXON

This model is described in detail in:

McIntyre CC, Richardson AG, and Grill WM. Modeling the excitability of
mammalian nerve fibers: influence of afterpotentials on the recovery
cycle. Journal of Neurophysiology 87:995-1006, 2002.
----------------------------------------------------------------------*/'''

def createMyelinatedParaFits():

    # we will extrapolate the McIntyre data in order to obtain arbitrary myelinated axon diameters

    # these values come directly from the McIntyre publication
    diams = [5.7, 7.3, 8.7, 10., 11.5, 12.8, 14., 15., 16.]
    nodeNodeSep = [500, 760, 1000, 1150, 1250, 1350, 1400, 1450, 1500]
    noLamella = [80, 100, 110, 120, 130, 135, 140, 145, 150]
    nodeDiam = [1.9, 2.4, 2.8, 3.3, 3.7, 4.2, 4.7, 5., 5.5]
    # MYSADiam = nodeDiam
    FLUTLen = [35, 38, 40, 46, 50, 54, 56, 58, 60]
    FLUTDiam = [3.4, 4.6, 5.8, 6.9, 8.1, 9.2, 10.4, 11.5, 12.7]
    STINLen = [70.5, 111.2, 152.2, 175.2, 190.5, 205.8, 213.5, 221.2, 228.8]
    # STIDiam = FLUTDiam

    # diamsMYSA = [0., 1., 1., 1., 16., 30, 40, 60, 80, 100.]
    # MYSALen = np.multiply([0.4, 0.7, 0.7, 0.7, 1, 1, 1, 1, 1, 1], 2.9)
    #
    # zMYSALen = np.polyfit(diamsMYSA, MYSALen, 5)
    # pMYSALen = np.poly1d(zMYSALen)

    shrinkingFactorMYSA = 0.4
    shrinkingFactorNode = shrinkingFactorMYSA

    diamsMYSA = [0., 16.]
    MYSALen = np.multiply([shrinkingFactorMYSA, 1], 3)

    zMYSALen = np.polyfit(diamsMYSA, MYSALen, 1)
    pMYSALen = np.poly1d(zMYSALen)

    diamsNode = [0., 16.]
    NodeLen = np.multiply([shrinkingFactorNode, 1], 1)

    # # new
    #
    # diamsMYSA = [           0.,     0.25,   0.5,    0.75,   1.,     16.]
    # MYSALen = np.multiply([ 0.15,    0.55,   0.6,    0.7,    0.8,    1], 3)
    #
    # MYSASpline = interpolate.splrep(diamsMYSA, MYSALen, k=1)
    #
    # diamsNode =             [0.,    0.25,   0.5,    0.75,   1.,     16.]
    # NodeLen = np.multiply(  [0.15,   0.45,   0.55,   0.6,    0.7,    1], 1)
    #
    # NodeSpline = interpolate.splrep(diamsNode, NodeLen, k=1)


    zNodeLen = np.polyfit(diamsNode, NodeLen, 1)
    pNodeLen = np.poly1d(zNodeLen)

    zNodeSep = np.polyfit(diams, nodeNodeSep, 1)
    pNodeSep = np.poly1d(zNodeSep)

    zNoLamella = np.polyfit(diams, noLamella, 1)
    pNoLamella = np.poly1d(zNoLamella)

    # diameters are fitted quadratically because linear fit results in negative values
    zNodeDiam = np.polyfit(diams, nodeDiam, 2)
    pNodeDiam = np.poly1d(zNodeDiam)

    zFLUTLen = np.polyfit(diams, FLUTLen, 1)
    pFLUTLen = np.poly1d(zFLUTLen)

    zSTINLen = np.polyfit(diams, STINLen, 1)
    pSTINLen = np.poly1d(zSTINLen)

    # diameters are fitted quadratically because linear fit results in negative values
    zFLUTDiam = np.polyfit(diams, FLUTDiam, 2)
    pFLUTDiam = np.poly1d(zFLUTDiam)

    return pNodeSep, pNoLamella, pNodeDiam, pFLUTLen, pSTINLen, pFLUTDiam, pMYSALen, pNodeLen # , MYSASpline, NodeSpline


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

    # static variables containing the fitted functions for diameter, node distances, etc.
    pNodeSep, pNoLamella, pNodeDiam, pFLUTLen, pSTINLen, pFLUTDiam, pMYSALen, pNodeLen = createMyelinatedParaFits() # , MYSASpline, nodeSpline

    @staticmethod
    def getFittedMcIntyreParams(diameter):

        axonD=Myelinated.pFLUTDiam(diameter)
        nodeD=Myelinated.pNodeDiam(diameter)
        paraD1=nodeD
        paraD2=axonD
        deltax=round(Myelinated.pNodeSep(diameter))
        paralength2=Myelinated.pFLUTLen(diameter)
        nl=np.round(Myelinated.pNoLamella(diameter))
        # interlength=int(round(Myelinated.pSTINLen(diameter)))
        g=axonD/diameter

        return axonD, nodeD, paraD1, paraD2, deltax, paralength2, nl,  g # interlength,

    @staticmethod
    def getOriginalMcIntyreParams(diameter):

        # find diameter closest to give nones in McIntyre's paper
        possibleDiams = np.array([5.7, 7.3, 8.7, 10., 11.5, 12.8, 14., 15., 16.])
        diffArray = np.abs(np.subtract(possibleDiams, diameter))
        diameterIndex = diffArray.argmin()

        fiberD = possibleDiams[diameterIndex]

        if fiberD==5.7:
            g=0.605
            axonD=3.4
            nodeD=1.9
            paraD1=1.9
            paraD2=3.4
            deltax=500
            paralength2=35
            nl=80
        elif fiberD==7.3:
            g=0.630
            axonD=4.6
            nodeD=2.4
            paraD1=2.4
            paraD2=4.6
            deltax=750
            paralength2=38
            nl=100
        elif fiberD==8.7:
            g=0.661
            axonD=5.8
            nodeD=2.8
            paraD1=2.8
            paraD2=5.8
            deltax=1000
            paralength2=40
            nl=110
        elif fiberD==10.0:
            g=0.690
            axonD=6.9
            nodeD=3.3
            paraD1=3.3
            paraD2=6.9
            deltax=1150
            paralength2=46
            nl=120
        elif fiberD==11.5:
            g=0.700
            axonD=8.1
            nodeD=3.7
            paraD1=3.7
            paraD2=8.1
            deltax=1250
            paralength2=50
            nl=130
        elif fiberD==12.8:
            g=0.719
            axonD=9.2
            nodeD=4.2
            paraD1=4.2
            paraD2=9.2
            deltax=1350
            paralength2=54
            nl=135
        elif fiberD==14.0:
            g=0.739
            axonD=10.4
            nodeD=4.7
            paraD1=4.7
            paraD2=10.4
            deltax=1400
            paralength2=56
            nl=140
        elif fiberD==15.0:
            g=0.767
            axonD=11.5
            nodeD=5.0
            paraD1=5.0
            paraD2=11.5
            deltax=1450
            paralength2=58
            nl=145
        else: #  (fiberD==16.0):
            g=0.791
            axonD=12.7
            nodeD=5.5
            paraD1=5.5
            paraD2=12.7
            deltax=1500
            paralength2=60
            nl=150

        return axonD, nodeD, paraD1, paraD2, deltax, paralength2, nl, g

    def __init__(self, fiberD, coord, tStop, timeRes, numberOfSavedSegments, temperature=37, v_init=-81, name="myelinated_axonA", layout3D="PT3D", rec_v=True, gkbar_axnode=0.12):
        super(Myelinated,self).__init__(layout3D, rec_v, name, fiberD, coord, temperature, v_init, tStop, timeRes, numberOfSavedSegments)

        self.endOverlap = 7 # number of axon segments after last node

        print 'Myelinated fiber diameter: ' + str(self.fiberD)

        # initialize variables to store NEURON elements in
        self.nodes = None
        self.MYSAs = None
        self.FLUTs = None
        self.STINs = None

        # morphological parameters
        self.paralength1 = 3 # MYSA length
        self.nodelength = 1.0 # node length
        self.space_p1 = 0.002
        self.space_p2 = 0.004
        self.space_i = 0.004

        #electrical parameters
        self.rhoa = 0.7e6 # Ohm-um
        self.mycm = 0.1 # uF/cm2/lamella membrane
        self.mygm = 0.001 # S/cm2/lamella membrane
        self.gkbar_axnode = gkbar_axnode

        # inter-/extrapolate the parameters of the McIntyre (2002) model
        axonD, nodeD, paraD1, paraD2, deltax, paralength2, nl, g = self.getFittedMcIntyreParams(self.fiberD) # interlength,

        # geometry including lammellae
        self.g = g
        self.axonD = axonD
        self.nodeD = nodeD
        self.paraD1 = paraD1
        self.paraD2 = paraD2
        self.deltax = deltax
        self.paralength2 = paralength2
        self.nl = nl

        # some additional length measures
        self.lengthOneCycle = self.deltax
        self.interlength = (self.deltax - self.nodelength - (2 * self.paralength1) - (2 * self.paralength2)) / 6

        # length of the whole axon
        self.L = createGeometry.length_from_coords(coord)

        # resistances
        self.Rpn0 = (self.rhoa*0.01)/(math.pi*((math.pow((self.nodeD/2)+self.space_p1,2))-(math.pow(nodeD/2,2))))
        self.Rpn1 = (self.rhoa*0.01)/(math.pi*((math.pow((self.paraD1/2)+self.space_p1,2))-(math.pow(paraD1/2,2))))
        self.Rpn2 = (self.rhoa*0.01)/(math.pi*((math.pow((self.paraD2/2)+self.space_p2,2))-(math.pow(paraD2/2,2))))
        self.Rpx = (self.rhoa*0.01)/(math.pi*((math.pow((self.axonD/2)+self.space_i,2))-(math.pow(axonD/2,2))))

        # number of sections
        # self.axonnodes = int(math.ceil(self.L/self.lengthOneCycle))
        self.axonnodes = int(math.ceil(self.L/self.lengthOneCycle)) - 1 # substract one for overlap
        self.paranodes1= 2*(self.axonnodes-1)
        self.paranodes2= 2*(self.axonnodes-1)
        self.axoninter= 6*(self.axonnodes-1)
        # self.axontotal= self.axonnodes+self.paranodes1+self.paranodes2+self.axoninter
        self.axontotal= self.axonnodes+self.paranodes1+self.paranodes2+self.axoninter + self.endOverlap

        print 'Number of segments for myelinated axon: %i' % self.axontotal

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

            # Attention, this is experimental. But e.g. in Roper and Schwarz 1989 'paranode regions' should also comprise MYSA segments in this model.
            # MYSA.insert('axflut') # self-made potassium channel with parameters from McIntyre et al. 2014 CNS data

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

            FLUT.insert('axflut') # self-made potassium channel with parameters from McIntyre et al. 2014 CNS data

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

            # # CAREFUL, this does not really belong here.
            # STIN.insert('axflut') # self-made potassium channel with parameters from McIntyre et al. 2014 CNS data

            self.STINs.append(STIN)
            self.allseclist.append(sec=STIN)

    def set_segment_types(self):

        nodeSequence = ['n', 'm', 'f', 's', 's', 's', 's', 's', 's', 'f', 'm'] # better be a member variable

        self.segmentTypes = []

        segCounter = 0
        for sec in self.allseclist:
            segmentTypeStr = nodeSequence[np.mod(segCounter, len(nodeSequence))]
            for seg in sec:
                self.segmentTypes.append(segmentTypeStr)
            segCounter += 1

    def create_neuron_object(self):

        self.nodes = []
        self.MYSAs = []
        self.FLUTs = []
        self.STINs = []

        # iterate through nodes in the way they are ordered in the axon in order to have an ordered allsectionlist
        nodeSequence = ['n', 'm', 'f', 's', 's', 's', 's', 's', 's', 'f', 'm']
        endSequence = nodeSequence[0:self.endOverlap+1] # add some
        for i in range(self.axonnodes-1):
            for j in range(len(nodeSequence)):
                nodeType = nodeSequence[j]
                self.createSingleNode(nodeType)
        for j in range(len(endSequence)):
            nodeType = endSequence[j]
            self.createSingleNode(nodeType)

        # transfer the geometry into the NEURON axon model
        self.position_sections_in_neuron()

        # do everything (number of segs per sec, channels, geometry from NEURON)
        super(Myelinated, self).create_neuron_object()

        # as number of segments per section are set dynamically depending on diameter, set the association here for
        # further processing
        self.set_segment_types()

        # here we correct the conductance of the slow potassium channel from 0.08 S/cm2 to 0.12 S/cm2 to prevent
        # multiple action potentials for thin fibers
        hString = 'forall for (x,0) if (ismembrane("axnode")) gkbar_axnode(x) = %2.3f' % self.gkbar_axnode
        h(hString)

    def delete_neuron_object(self):

        super(Myelinated, self).delete_neuron_object()

        self.nodes = None
        self.FLUTs = None
        self.MYSAs = None
        self.STINs = None