import numpy as np
import time
from PyPN.takeTime import *
from scipy.interpolate import Rbf
from scipy.interpolate import NearestNDInterpolator
import scipy.ndimage
import fractions
from scipy import signal
import sys
import matplotlib.pyplot as plt

def show_sizeof(x, level=0):

    print "\t" * level, x.__class__, sys.getsizeof(x), x

    if hasattr(x, '__iter__'):
        if hasattr(x, 'items'):
            for xx in x.items():
                show_sizeof(xx, level + 1)
        else:
            for xx in x:
                show_sizeof(xx, level + 1)

def downsample(s, n, phase=0):
    """Decrease sampling rate by integer factor n with included offset phase.
    """
    if len(np.shape(s)) == 1:
        return s[phase::n]
    else:
        return s[:,phase::n]

def upsample(s, n, phase=0):
    """Increase sampling rate by integer factor n  with included offset phase.
    """
    return np.roll(np.kron(s, np.r_[1, np.zeros(n-1)]), phase)

def upfirdn(s, h, p, q):
    """Upsample signal s by p, apply FIR filter as specified by h, and
    downsample by q. Using fftconvolve as opposed to lfilter as it does not seem
    to do a full convolution operation (and its much faster than convolve).
    """

    us = upsample(s, p)
    lowpassed = scipy.ndimage.convolve1d(us, h, axis=-1)
    return downsample(lowpassed, q)
    # return downsample(signal.fftconvolve(np.tile(h, np.shape(us)), us), q)
    # return downsample(upsample(s, p), q)

def resample(s, p, q, h=None):
    """Change sampling rate by rational factor. This implementation is based on
    the Octave implementation of the resample function. It designs the
    anti-aliasing filter using the window approach applying a Kaiser window with
    the beta term calculated as specified by [2].

    Ref [1] J. G. Proakis and D. G. Manolakis,
    Digital Signal Processing: Principles, Algorithms, and Applications,
    4th ed., Prentice Hall, 2007. Chap. 6

    Ref [2] A. V. Oppenheim, R. W. Schafer and J. R. Buck,
    Discrete-time signal processing, Signal processing series,
    Prentice-Hall, 1999
    """
    gcd = fractions.gcd(p, q)
    if gcd > 1:
        p = p / gcd
        q = q / gcd

    if h is None:  # design filter
        # properties of the antialiasing filter
        log10_rejection = -3.0
        stopband_cutoff_f = 1.0 / (2.0 * max(p, q))
        roll_off_width = stopband_cutoff_f / 10.0

        # determine filter length
        # use empirical formula from [2] Chap 7, Eq. (7.63) p 476
        rejection_db = -20.0 * log10_rejection;
        l = np.ceil((rejection_db - 8.0) / (28.714 * roll_off_width))

        # ideal sinc filter
        t = np.arange(-l, l + 1)
        ideal_filter = 2 * p * stopband_cutoff_f * np.sinc(2 * stopband_cutoff_f * t)

        # determine parameter of Kaiser window
        # use empirical formula from [2] Chap 7, Eq. (7.62) p 474
        beta = signal.kaiser_beta(rejection_db)

        # apodize ideal filter response
        h = np.kaiser(2 * l + 1, beta) * ideal_filter

    ls = np.shape(s)[-1]
    lh = len(h)

    l = (lh - 1) / 2.0
    ly = np.ceil(ls * p / float(q))

    # pre and postpad filter response
    nz_pre = np.floor(q - np.mod(l, q))
    hpad = h[-lh + nz_pre:]

    offset = np.floor((l + nz_pre) / q)
    nz_post = 0
    while np.ceil(((ls - 1) * p + nz_pre + lh + nz_post) / q) - offset < ly:
        nz_post += 1
    hpad = hpad[:lh + nz_pre + nz_post]

    # filtering
    xfilt = upfirdn(s, hpad, p, q)

    if len(np.shape(xfilt)) == 1:
        return xfilt[0:ly]
    elif len(np.shape(xfilt)) == 2:
        return xfilt[:, 0:ly]
    else:
        return []

    # return xfilt[offset - 1:offset - 1 + ly]

def rectangularStimulusSignal(stimDur, amplitude, frequency, dutyCycle, waveform, timeRes, delay=0, invert=False):

    t = np.arange(0, stimDur, timeRes)

    if waveform == 'MONOPHASIC':
        stimulusSignal = amplitude * 0.5 * signal.square(2 * np.pi * frequency * t, duty=dutyCycle) + amplitude * 0.5
    elif waveform == 'BIPHASIC':
        stimulusSignal = amplitude * signal.square(2 * np.pi * frequency * t, duty=dutyCycle)
    else:
        print "You didn't choose the right waveform either MONOPHASIC or BIPHASIC, it has been set to default MONOPHASIC"
        stimulusSignal = amplitude * 0.5 * signal.square(2 * np.pi * frequency * t, duty=dutyCycle) + amplitude * 0.5

    if invert:
        stimulusSignal = -stimulusSignal

    return t, stimulusSignal

# t, signalRec = rectangularStimulusSignal(100, 1, 1, 0.01, 'MONOPHASIC', 0.5)

def change_samplingrate(inputSignal, inputTime, factor):

    # translate factor to rational ratio
    pq = fractions.Fraction(str(factor)).limit_denominator(max_denominator=100)
    p = pq.numerator
    q = pq.denominator

    # resampled signal
    signalResampled = resample(inputSignal, p, q)

    # resampled time
    dt = np.mean(np.diff(inputTime))
    dtResampled = dt/factor
    tResampled = np.arange(0, np.max(inputTime), dtResampled)
    lendiff = len(tResampled) - np.shape(signalResampled)[-1]
    if lendiff > 0:
        len(tResampled) > len(signalResampled)
        tResampled = tResampled[:-(lendiff)]
    elif lendiff < 0:
        tResampled = np.concatenate((tResampled, tResampled[-1] + (dtResampled * np.arange(1, -lendiff + 1))))

    return tResampled, signalResampled

# # resampling factor
# factor = 31.0
# pq = fractions.Fraction(str(factor)).limit_denominator(max_denominator=100)
# p = pq.numerator
# q = pq.denominator
#
# # sample signal
# t = np.arange(0,1000,.1)
# signalRec = np.sin(2*np.pi/1000*t)
#
# # resampled signal
# signalRecResampled = resample(signalRec, p, q)
#
# # resampled time
# dt = np.mean(np.diff(t))
# dtResampled = dt/factor
# tResampled = np.arange(0,max(t), dtResampled)
# lendiff = len(tResampled) - len(signalRecResampled)
# if lendiff > 0:
#     len(tResampled) > len(signalRecResampled)
#     tResampled = tResampled[:-(1+lendiff)]
# elif lendiff < 0:
#     tResampled = np.concatenate((tResampled, tResampled[-1]+(dtResampled*np.arange(1, -lendiff+1))))

# # sample signal
# t = np.arange(0,1000,.1)
# signalRec = np.sin(2*np.pi/1000*t)

t = np.tile(np.arange(0,1000,.1), (2, 1))
signalRec = np.sin(2*np.pi/1000*t)

tResampled, signalRecResampled = change_samplingrate(signalRec, t, 0.05)

plt.plot(t[0,:], signalRec[0,:], '-', color=(0,1,0))#, linestyle='steps--')
# plt.plot(tResampled, signalRecResampled[0,:], '--', color=(1,0,0))#, linestyle='steps--') # tResampled,
plt.plot(tResampled, signalRecResampled[1,:], '--', color=(1,0,0))#, linestyle='steps--') # tResampled,
# plt.stem(t, signalRec, color=(0,1,0))
plt.show()


