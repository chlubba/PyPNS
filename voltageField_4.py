"""
To interpolate voltage,
1. get entire field for every time step
2. interpolate for all points of interest (normally segment positions) and append them to vectors

Always keep track of number of points


"""

import numpy as np
import time
from PyPN.takeTime import *
from scipy.interpolate import Rbf, griddata, NearestNDInterpolator, LinearNDInterpolator
import fractions
from scipy import signal
import sys

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
    return s[phase::n]


def upsample(s, n, phase=0):
    """Increase sampling rate by integer factor n  with included offset phase.
    """
    return np.roll(np.kron(s, np.r_[1, np.zeros(n-1)]), phase)

def upfirdn(s, h, p, q):
    """Upsample signal s by p, apply FIR filter as specified by h, and
    downsample by q. Using fftconvolve as opposed to lfilter as it does not seem
    to do a full convolution operation (and its much faster than convolve).
    """
    return downsample(signal.fftconvolve(h, upsample(s, p)), q)

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

    ls = len(s)
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

    return xfilt[offset - 1:offset - 1 + ly]

def printSize(array):
    nbytes = array.nbytes
    mbytes = np.floor(nbytes / (1024 ** 2))
    remainder = nbytes - 1024 ** 2 * mbytes
    kbytes = np.floor(remainder / 1024)
    bytes = remainder - 1024 * kbytes
    print '%i MB %i KB %i B' % (mbytes, kbytes, bytes)

# ------------- body starts ---------------

fieldPoints = 100000
segmentPoints = 10000
timeSteps = 1

# check if

xi = yi = zi = np.linspace(0, 1, segmentPoints)
pointsOfInterest = np.vstack([xi, yi, zi])

interpolatedVoltageTimeTemp = np.array([]).reshape(0,pointsOfInterest.shape[1])

t0 = time.time()
for i in range(timeSteps):

    voltageCoordsTemp = np.random.rand(fieldPoints, 3)
    voltageValuesTemp = np.squeeze(np.random.rand(1, fieldPoints))

    # x_ = y_ = z_ = np.linspace(0, 2000, 100)
    #
    # xx, yy, zz = np.meshgrid(x_, y_, z_)
    #
    # x = np.squeeze(xx.reshape(1, -1))
    # y = np.squeeze(yy.reshape(1, -1))
    # z = np.squeeze(zz.reshape(1, -1))
    #
    # voltageValuesTemp = np.sin(2 * np.pi / 2000 * x) * np.sin(2 * np.pi / 2000 * y) * np.sin(2 * np.pi / 2000 * y)
    # voltageCoordsTemp = np.vstack([x, y, z])


    printSize(voltageCoordsTemp)

    with takeTime('interpolator'):

        interpolatorN = NearestNDInterpolator(voltageCoordsTemp, voltageValuesTemp)
        interpoltorLin = LinearNDInterpolator(voltageCoordsTemp, voltageValuesTemp)
    # interpolatorRbf = Rbf(voltageCoordsTemp[:,0], voltageCoordsTemp[:,1], voltageCoordsTemp[:,2], voltageValuesTemp)  # radial basis function interpolator instance

    interpolatedVoltageTemp = []
    # with takeTime('interpolation'):
    # go through all segment positions for the current time step and save values
    # interpolatedValues = interpolatorN(pointsOfInterest[0,:], pointsOfInterest[1,:], pointsOfInterest[2,:])  # interpolated values
    for poiIndex in range(pointsOfInterest.shape[1]):
        poi = np.squeeze(pointsOfInterest[:,poiIndex])
        diN = interpolatorN(poi[0], poi[1], poi[2])   # interpolated values
        # diRbf = interpolatorRbf(xi, yi, zi)  # interpolated values
        interpolatedVoltageTemp.append(diN)

    interpolatedVoltageTimeTemp = np.vstack([interpolatedVoltageTimeTemp, interpolatedVoltageTemp])

    # nbytes = interpolatedVoltageTimeTemp.nbytes
    # mbytes = np.floor(nbytes/(1024**2))
    # remainder = nbytes - 1024**2*mbytes
    # kbytes = np.floor(remainder/1024)
    # bytes = remainder - 1024 * kbytes
    #
    if i == 0:
        # print 'Per time step: %i MB %i KB %i B' % (mbytes, kbytes, bytes)
        print 'One time step: ',
        printSize(interpolatedVoltageTimeTemp)

# nbytes = interpolatedVoltageTimeTemp.nbytes
# mbytes = np.floor(nbytes/(1024**2))
# remainder = nbytes - 1024**2*mbytes
# kbytes = np.floor(remainder/1024)
# bytes = remainder - 1024 * kbytes
# print 'Overall: %i MB %i KB %i B' % (mbytes, kbytes, bytes)
print 'Overall size of interpolated values on disk:',
printSize(interpolatedVoltageTimeTemp)

print 'Overall time for %i field points, %i time steps, %i points of interest: %3.2f' \
      % (fieldPoints, timeSteps, segmentPoints, time.time() - t0)


