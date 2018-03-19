from neuron import h
import numpy as np
import math
from scipy import signal
from excitationMechanismClass import *
import os

from signalGeneration import *
from samplingRates import *
from nameSetters import *
from extracellularBackend import *


class StimIntra(ExcitationMechanism):
    # stimType: INTRA or EXTRA cellular stimulation
    # axon: axon object on which the stimulation is applied
    # pos: position of the stimulus
    # sect: section being stimulated
    # delay: pulse delay (ms)
    # stimDur: pulse duration (ms)
    # amp: pulse amplitude (nA)
    # freq: frequency of the sin pulse (Hz)
    # duty_cycle: Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
    # stim_coord=[xe,ye,ze]: spatial coordinates  of the stimulating electrode
    # waveform: Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric

    def __init__(self, stimulusSignal):
        """Intracellular stimulation is effectively a current injection into each axon.

        :param stimulusSignal: signal to stimulate the axons with. Unit nA
        """

        self.stimulusSignal = stimulusSignal
        self.svec = h.Vector(self.stimulusSignal)

        super(StimIntra, self).__init__()

    def connect_axon(self, axon):
        """ An IClamp is placed on the first section of the axon.

        :param axon: The axon the mechanism is connected to.

        """

        # Place an IClamp on the first element of the allseclist
        # In unmyelinated axon case allseclist is directly the unique axon section

        stim = h.IClamp(0, axon.allseclist)
        stim.dur = len(self.stimulusSignal)*self.timeRes
        self.svec.play(stim._ref_amp, self.timeRes)

        excitationMechanismVars = [stim]
        axon.append_ex_mech_vars(excitationMechanismVars)

    def delete_neuron_objects(self):
        self.svec = None


class StimField(ExcitationMechanism):
    def __init__(self, stimulusSignal, electrodePositions, extPotMech, polarities = ()):
        """A field is calulated from the current on stimulation electrodes and then applied to the axon as the extracellular membrane voltage.

        :param stimulusSignal: Current signal on the electrodes [nA]. If electrode is composed of many point electrodes, the signal is divided by the number of those.
        :param electrodePositions: Positions of the electrodes as for RecordingMechanisms 3D
        :param extPotMech: ``ExtracellularMechansim`` used to calculate the voltage at the axon membrane caused by the stimulation current.
        :param polarities: Signs of the stimulation electrodes.
        """
        super(StimField, self).__init__()

        # electrode setup: positions and polarities
        self.electrodePositions = electrodePositions
        self.numberOfPoints = np.shape(electrodePositions)[0]
        self.numberOfPoles = np.shape(electrodePositions)[2]
        if polarities == ():
            self.polarities = np.power((-1), range(self.numberOfPoles))
        else:
            self.polarities = polarities
        assert (self.numberOfPoles == len(self.polarities))

        self.extPotMech = extPotMech

        self.signal = stimulusSignal

    def connect_axon(self, axon):

        # one signal for every point that constitutes an electrode. Divide by number of points to keep current constant
        signalTiled = np.tile(self.signal/self.numberOfPoints, (self.numberOfPoints, 1))

        # calculate the potential caused by stimulation poles for reference current of constant 1nA, sum everything up
        extSegPot = np.zeros((len(axon.xmid), 1))
        for poleIndex in range(self.numberOfPoles):
            # todo: check current units, and voltage units...
            extSegPot += self.polarities[poleIndex] * \
                         self.extPotMech.calculate_extracellular_potential(self.electrodePositions[:, :, poleIndex],
                                                                           np.ones(np.shape(self.electrodePositions)[0])[:, np.newaxis],
                                                                           np.transpose(np.vstack([axon.xmid, axon.ymid, axon.zmid])))

        # interpret calculated voltage as a transfer resistivity r = v_ext/i_ref with i_ref = 1nA
        segCounter = 0
        for sec in axon.allseclist:
            for segInd, seg in enumerate(sec):
                seg.xtra.rx = extSegPot[segCounter] * 1e-6 # see xtra readme online from neuron.yale.org
                # print 'seg no %i, rx %20.20f' % (segCounter, seg.xtra.rx)
                segCounter += 1

        # write signal vector into the reference for the xtra mechanism
        svec = h.Vector(np.concatenate((self.signal, [0])))
        svec.play(h._ref_is_xtra, constants.timeResStim)

        # keep variables in memory in order for NEURON to see them
        axon.append_ex_mech_vars([svec])


    def delete_neuron_objects(self):
        pass

class SimpleIClamp(ExcitationMechanism):

    def __init__(self, delay, stimDur, amplitude):

        self.delay = delay
        self.stimDur = stimDur
        self.amplitude = amplitude

        super(SimpleIClamp, self).__init__()


    def connect_axon(self, axon):

        # Place an IClamp on the first element of the allseclist
        # In unmyelinated axon case allseclist is directly the unique axon section

        stim = h.IClamp(0.001, axon.allseclist)
        stim.delay = self.delay
        stim.dur = self.stimDur
        stim.amp = self.amplitude

        excitationMechanismVars = [stim]
        axon.append_ex_mech_vars(excitationMechanismVars)

    def delete_neuron_objects(self):
        pass
