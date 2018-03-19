import numpy as np
from scipy.interpolate import interp1d

def idealizedCuff(cuffWidthIn):

    cuffWidth = cuffWidthIn/2

    a = 2.5E-9 # 1.9E-9  # delta peak max
    b = 0.00005
    nerveWidth = 0.00019
    cuffWidth0 = 0.01
    triangleMax0 = 8.83e-5

    # with decreased length the resistance also decreases and therefore the voltage
    sizeFactor = cuffWidth/cuffWidth0
    triangleMax = triangleMax0 * sizeFactor

    # for z-dependent triangle, use interpolation
    def smooth(y, box_pts):
        box = np.ones(box_pts) / box_pts
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth

    dz = min((cuffWidth/500,0.0001))
    interpLen = 1.5*cuffWidth
    zInterp = np.arange(-interpLen, interpLen, dz)
    if len(zInterp)%2==1:
        zInterp=zInterp[0:-1]
    smoothWidth = cuffWidth/20 # 5
    smoothSamples = smoothWidth / dz
    sharpOneSide = np.maximum(0, triangleMax * (np.add(1, np.divide(zInterp, cuffWidth))))
    smoothedOneSide = smooth(sharpOneSide, int(smoothSamples))
    smoothedOneSideToMiddle = smoothedOneSide[0:int(np.floor(np.shape(smoothedOneSide)[0] / 2))]
    smoothedTwoSides = np.concatenate([smoothedOneSideToMiddle, np.fliplr([smoothedOneSideToMiddle])[0]])
    triangle = interp1d(zInterp, smoothedTwoSides, bounds_error=False, fill_value=0)

    peakFactor = lambda angle, xP: np.maximum(0, (1 - np.abs(np.mod(angle + np.pi, 2*np.pi)-np.pi) / np.pi * 5)) * np.minimum(
        1, (xP / nerveWidth) ** 5)
    peak = lambda zValues, angle, xP: a * (1.0 / (np.abs(zValues) + b)) * peakFactor(angle, xP)

    return lambda zValues, angle, xP: triangle(zValues) + peak(zValues, angle, xP)