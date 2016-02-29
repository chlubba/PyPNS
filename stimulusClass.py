from neuron import h
import numpy as np # for arrays managing
import math

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