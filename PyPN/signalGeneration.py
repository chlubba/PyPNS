import numpy as np
from scipy import signal
import constants

def biphasic_decaying(timeRes=constants.timeResStim, tDelay=0, tC=1, aC=0.05, tExp=1, cExp=-5, tD=2, aD=-0.01):
    """
    Defines stimulus signal shape as in <name paper>

    Args:
        timeRes: time resolution of the simulation
        tDelay: length of first zero phase [ms]
        tC: length of up phase [ms]
        aC: amplitude of up phase [mA]
        tExp: length of decaying phase
        cExp: decaying coefficient
        tD: length down time
        aD: amplitude down time

    Returns:
        t: time array
        overallSignal: signal array

    """


    signal0 = np.zeros(int(tDelay / timeRes))
    signal1 = np.ones(int(tC / timeRes)) * aC
    tTempExp = np.arange(0, tExp, timeRes)
    signal2 = np.exp(tTempExp * cExp) * aC
    signal3 = np.ones(int(tD / timeRes)) * aD

    overallSignal = np.concatenate((signal0, signal1, signal2, signal3))

    # finalize signal with a zero
    overallSignal = np.concatenate((overallSignal, [0]))

    # t = np.arange(len(overallSignal)) * tRes

    return overallSignal

def rectangular(stimDur, amplitude, frequency, dutyCycle, waveform, timeRes=constants.timeResStim, delay=0, invert=False):

    tGen = np.arange(0, stimDur, timeRes)

    if waveform == 'MONOPHASIC':
        stimulusSignal = amplitude * 0.5 * signal.square(2 * np.pi * frequency * tGen, duty=dutyCycle) + amplitude * 0.5
    elif waveform == 'BIPHASIC':
        stimulusSignal = amplitude * signal.square(2 * np.pi * frequency * tGen, duty=dutyCycle)
    else:
        print "You didn't choose the right waveform either MONOPHASIC or BIPHASIC, it has been set to default MONOPHASIC"
        stimulusSignal = amplitude * 0.5 * signal.square(2 * np.pi * frequency * tGen, duty=dutyCycle) + amplitude * 0.5

    if invert:
        stimulusSignal = -stimulusSignal

    # finalize signal with a zero
    stimulusSignal = np.concatenate((stimulusSignal, [0]))

    # add delay
    stimulusSignalDelayed = np.concatenate((np.zeros(delay / timeRes), stimulusSignal))
    # t = np.arange(len(stimulusSignalDelayed))*timeRes

    return stimulusSignalDelayed

