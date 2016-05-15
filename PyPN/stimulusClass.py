from neuron import h
import numpy as np
import math
from scipy import signal
from excitationMechanismClass import *

class Stimulus(ExcitationMechanism):
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
    def __init__(self, stimType, stimDur, amplitude, frequency, dutyCycle, radiusBundle, waveform, timeRes, delay=0):

        self.waveform = waveform
        self.stimType = stimType

        number_contact_points = 8
        self.stim_coord = [[0, radiusBundle*math.cos(math.pi/number_contact_points),
                            radiusBundle*math.sin(math.pi/number_contact_points)]]

        self.stimDur = stimDur
        self.frequency = frequency
        self.amplitude = amplitude
        self.dutyCycle = dutyCycle
        self.delay = delay

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

        # # add delay
        # self.stimulusSignalDelayed = np.concatenate((np.zeros(self.delay/timeRes), self.stimulusSignal))

        # cut_off = math.sin((-self.dutyCycle+1.0/2)*math.pi)
        # if self.waveform == "MONOPHASIC":
        #     self.signal = self.amplitude*(cut_off < np.sin(2*math.pi*self.frequency*(self.t)))
        # elif self.waveform == "BIPHASIC":
        #     self.signal = self.amplitude*(cut_off < np.sin(2*math.pi*self.frequency*self.t))-self.amplitude*(cut_off < np.sin(2*math.pi*self.frequency*(self.t+self.dutyCycle/self.frequency)))
        # else:
        #     print "You didn't choose the right waveform either MONOPHASIC or BIPHASIC, it has been set to default MONOPHASIC"
        #     self.signal = self.amplitude*(cut_off < np.sin(2*math.pi*self.frequency*(self.t)))

        self.svec = h.Vector(self.stimulusSignal)

    def connect_axon(self, axon):
        if self.stimType == "INTRA":
            self.init_intra(axon)
        elif self.stimType == "EXTRA":
            axon.setrx(self.stim_coord, axon.axonPosition)
            self.init_xtra()
        else:
            raise NameError('stimType only "INTRA" or "EXTRA"')

    def init_intra(self, axon):
        # Place an IClamp on the first element of the allseclist
        # In unmyelinated axon case allseclist is directly the unique axon section

        stim = h.IClamp(0, axon.allseclist)
        stim.delay = self.delay
        stim.dur = self.stimDur
        self.svec.play(stim._ref_amp, self.timeRes)

        excitationMechanismVars = [stim]
        axon.append_ex_mech_vars(excitationMechanismVars)

    def init_xtra(self):
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