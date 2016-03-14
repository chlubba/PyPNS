from neuron import h
import numpy as np
import math
from ExcitationMechanism import *

class Stimulus(ExcitationMechanism):
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
    def __init__(self, stim_type, dur, amp, freq, duty_cycle, radius_bundle, waveform):
        self.waveform = waveform
        self.stim_type = stim_type

        number_contact_points = 8
        self.stim_coord = [[0, radius_bundle*math.cos(math.pi/number_contact_points),
                            radius_bundle*math.sin(math.pi/number_contact_points)]]
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

    def connect_axon(self, axon):
        if self.stim_type == "INTRA":
            self.init_intra(axon)
        elif self.stim_type == "EXTRA":
            axon.setrx(self.stim_coord, axon.axonPosition)
            self.init_xtra()
        else:
            raise NameError('stim_type only "INTRA" or "EXTRA"')

    def init_intra(self, axon):
        # Place an IClamp on the first element of the allseclist
        # In unmyelinated axon case allseclist is directly the unique axon section

        stim = h.IClamp(0, axon.allseclist)
        stim.delay = 0
        stim.dur = self.dur
        self.svec.play(stim._ref_amp,h.dt)

        excitationMechanismVars = [stim]
        axon.appendExMechVars(excitationMechanismVars)

    def init_xtra(self):
        self.svec.play(h._ref_is_xtra,h.dt)

    def delete_neuron_objects(self):
        self.svec = None