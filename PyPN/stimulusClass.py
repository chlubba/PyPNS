from neuron import h
import numpy as np
import math
from scipy import signal
from excitationMechanismClass import *



class StimIntra(ExcitationMechanism):
    """
    stimType: INTRA or EXTRA cellular stimulation
    axon: axon object on which the stimulation is applied
    pos: position of the stimulus
    sect: section being stimulated
    delay: pulse delay (ms)
    stimDur: pulse duration (ms)
    amp: pulse amplitude (nA)
    freq: frequency of the sin pulse (Hz)
    duty_cycle: Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
    stim_coord=[xe,ye,ze]: spatial coordinates  of the stimulating electrode
    waveform: Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
    """
    def __init__(self, stimDur, amplitude, frequency, dutyCycle, waveform, timeRes, delay=0, invert=False):

        self.waveform = waveform

        self.stimDur = stimDur
        self.frequency = frequency
        self.amplitude = amplitude
        self.dutyCycle = dutyCycle
        self.delay = delay

        self.timeRes = timeRes

        self.t = np.arange(0, self.stimDur, timeRes)

        if self.waveform == 'MONOPHASIC':
            self.stimulusSignal = self.amplitude * 0.5 * signal.square(2 * np.pi * self.frequency * self.t, duty=dutyCycle) + self.amplitude * 0.5
        elif self.waveform == 'BIPHASIC':
            self.stimulusSignal = self.amplitude * signal.square(2 * np.pi * self.frequency * self.t, duty=dutyCycle)
        else:
            print "You didn't choose the right waveform either MONOPHASIC or BIPHASIC, it has been set to default MONOPHASIC"
            self.stimulusSignal = self.amplitude * 0.5 * signal.square(2 * np.pi * self.frequency * self.t, duty=dutyCycle) + self.amplitude * 0.5

        if invert:
            self.stimulusSignal = -self.stimulusSignal


        # add delay
        self.stimulusSignalDelayed = np.concatenate((np.zeros(self.delay/timeRes), self.stimulusSignal))

        # end with a zero, otherwise final value is valid for the rest of the simulation...
        self.stimulusSignalDelayed = np.concatenate((self.stimulusSignalDelayed, [0]))

        self.svec = h.Vector(self.stimulusSignalDelayed)

    def connect_axon(self, axon):

        # Place an IClamp on the first element of the allseclist
        # In unmyelinated axon case allseclist is directly the unique axon section

        stim = h.IClamp(0, axon.allseclist)
        stim.dur = self.stimDur + self.delay
        self.svec.play(stim._ref_amp, self.timeRes)

        excitationMechanismVars = [stim]
        axon.append_ex_mech_vars(excitationMechanismVars)

    def delete_neuron_objects(self):
        self.svec = None

class StimCuff(ExcitationMechanism):
    """
    stimType: INTRA or EXTRA cellular stimulation
    axon: axon object on which the stimulation is applied
    pos: position of the stimulus
    sect: section being stimulated
    delay: pulse delay (ms)
    stimDur: pulse duration (ms)
    amp: pulse amplitude (nA)
    freq: frequency of the sin pulse (Hz)
    duty_cycle: Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
    stim_coord=[xe,ye,ze]: spatial coordinates  of the stimulating electrode
    waveform: Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
    """
    def __init__(self, stimDur, amplitude, frequency, dutyCycle, radius, waveform, timeRes, delay=0, invert=False):

        self.waveform = waveform

        numberContactPoints = 20 # more or less angle-independent field (<1% difference at 80% of the electrode radius)
        angles = 2*np.pi/numberContactPoints*np.arange(numberContactPoints)
        self.stim_coord = np.column_stack((np.zeros(numberContactPoints), radius*np.cos(angles), radius*np.sin(angles)))

        self.stimDur = stimDur
        self.frequency = frequency
        self.amplitude = amplitude/numberContactPoints
        self.dutyCycle = dutyCycle
        self.delay = delay
        self.radius = radius

        self.timeRes = timeRes

        # self.t = np.linspace(0, self.stimDur, (tStop+2*timeRes)/timeRes, endpoint=True)
        self.t = np.arange(0, self.stimDur, timeRes)

        if self.waveform == 'MONOPHASIC':
            self.stimulusSignal = self.amplitude * 0.5 * signal.square(2 * np.pi * self.frequency * self.t, duty=dutyCycle) + self.amplitude * 0.5
        elif self.waveform == 'BIPHASIC':
            self.stimulusSignal = self.amplitude * signal.square(2 * np.pi * self.frequency * self.t, duty=dutyCycle)
        else:
            print "You didn't choose the right waveform either MONOPHASIC or BIPHASIC, it has been set to default MONOPHASIC"
            self.stimulusSignal = self.amplitude * 0.5 * signal.square(2 * np.pi * self.frequency * self.t, duty=dutyCycle) + self.amplitude * 0.5

        if invert:
            self.stimulusSignal = -self.stimulusSignal


        # add delay
        self.stimulusSignalDelayed = np.concatenate((np.zeros(self.delay/timeRes), self.stimulusSignal))

        # end with a zero, otherwise final value is valid for the rest of the simulation...
        self.stimulusSignalDelayed = np.concatenate((self.stimulusSignalDelayed, [0]))

        self.svec = h.Vector(self.stimulusSignalDelayed)

    def connect_axon(self, axon):

        axon.setrx(self.stim_coord, axon.axonPosition)
        self.svec.play(h._ref_is_xtra, self.timeRes)


    def delete_neuron_objects(self):
        self.svec = None


class SimpleIClamp(ExcitationMechanism):

    def __init__(self, delay, stimDur, amplitude):

        self.delay = delay
        self.stimDur = stimDur
        self.amplitude = amplitude


    def connect_axon(self, axon):

        # Place an IClamp on the first element of the allseclist
        # In unmyelinated axon case allseclist is directly the unique axon section

        stim = h.IClamp(0, axon.allseclist)
        stim.delay = self.delay
        stim.dur = self.stimDur
        stim.amp = self.amplitude

        excitationMechanismVars = [stim]
        axon.append_ex_mech_vars(excitationMechanismVars)

    def delete_neuron_objects(self):
        pass