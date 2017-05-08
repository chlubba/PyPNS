import numpy as np
from scipy.interpolate import interp1d

def idealCuff(cuffWidth=0.01, r1=0.00019, a=2.5E-9, b=0.00005, triangleMax=8.83E-5):

    # a = 2.5E-9  # 1.9E-9  #
    # b = 0.00005
    # cuffWidth = 0.01
    # triangleMax = 8.83e-5

    # for z-dependent triangle, use interpolation
    def smooth(y, box_pts):
        box = np.ones(box_pts) / box_pts
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth

    dz = 0.000001
    zInterp = np.arange(-0.02, 0.02, dz)
    smoothWidth = cuffWidth / 5
    smoothSamples = smoothWidth / dz
    sharpOneSide = np.maximum(0, triangleMax * (np.add(1, np.divide(zInterp, cuffWidth))))
    smoothedOneSide = smooth(sharpOneSide, int(smoothSamples))
    smoothedOneSideToMiddle = smoothedOneSide[0:int(np.floor(np.shape(smoothedOneSide)[0] / 2))]
    smoothedTwoSides = np.concatenate([smoothedOneSideToMiddle, np.fliplr([smoothedOneSideToMiddle])[0]])
    triangle = interp1d(zInterp, smoothedTwoSides)

    peakFactor = lambda angle, xP: np.maximum(0, (
    1 - np.abs(np.mod(angle + np.pi, 2 * np.pi) - np.pi) / np.pi * 5)) * np.minimum(
        1, (xP / r1) ** 5)
    peak = lambda zValues, angle, xP: a * (1.0 / (np.abs(zValues) + b)) * peakFactor(angle, xP)

    return lambda zValues, angle, xP: triangle(zValues) + peak(zValues, angle, xP)