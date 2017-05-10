from abc import abstractmethod, ABCMeta
import numpy as np
import os
import silencer
import time
from extracellularBackend import *
from scipy.interpolate import interp1d

class ExtracellularPotentialMechanism(object):
    __metaclass__ = ABCMeta

    @abstractmethod

    def calculate_extracellular_potential(self, sourcePositions, sourceCurrents, receiverPositions):
        pass


class analytic(ExtracellularPotentialMechanism):

    def __init__(self, bundleGuide, method='z,xP,angle', interpolator=None):
        """Uses and idealized interpolation function to calculate the extracelluar potential caused by point sources. Used by :class:`recodingMechanismClass.RecordingMechanism`.

        :param bundleGuide: 3D nerve trajectory
        :param fieldName: string containing the name of the field to import; the location of field files needs to be specified somewhere (TODO)
        """

        self.bundleGuide = bundleGuide
        self.interpolator = interpolator
        self.method = method

        if self.interpolator == None:
            # if no function for potential calculation is given, take the default one

            a = 2.5E-9 # 1.9E-9  #
            b = 0.00005
            cuffWidth = 0.01
            triangleMax = 8.83e-5

            # for z-dependent triangle, use interpolation
            def smooth(y, box_pts):
                box = np.ones(box_pts) / box_pts
                y_smooth = np.convolve(y, box, mode='same')
                return y_smooth

            dz = 0.000001
            zInterp = np.arange(-0.02, 0.02, dz)
            smoothWidth = cuffWidth/20 # 5
            smoothSamples = smoothWidth / dz
            sharpOneSide = np.maximum(0, triangleMax * (np.add(1, np.divide(zInterp, cuffWidth))))
            smoothedOneSide = smooth(sharpOneSide, int(smoothSamples))
            smoothedOneSideToMiddle = smoothedOneSide[0:int(np.floor(np.shape(smoothedOneSide)[0] / 2))]
            smoothedTwoSides = np.concatenate([smoothedOneSideToMiddle, np.fliplr([smoothedOneSideToMiddle])[0]])
            triangle = interp1d(zInterp, smoothedTwoSides, bounds_error=False, fill_value="extrapolate")

            peakFactor = lambda angle, xP: np.maximum(0, (1 - np.abs(np.mod(angle + np.pi, 2*np.pi)-np.pi) / np.pi * 5)) * np.minimum(
                1, (xP / 0.000190) ** 5)
            peak = lambda zValues, angle, xP: a * (1.0 / (np.abs(zValues) + b)) * peakFactor(angle, xP)

            self.interpolator = lambda zValues, angle, xP: triangle(zValues) + peak(zValues, angle, xP)


            # import matplotlib.pyplot as plt
            # plt.figure()
            # # zPlot = np.arange(-0.01, 0.01, 0.0001)
            # angles = np.arange(-5*np.pi, 5*np.pi, 0.001)
            # plt.plot(angles, peakFactor(angles, 0.00019))
            # # angle = 0 # np.pi/10
            # # for xPPlot in [0.00009, 0.00018, 0.00025]:
            # #     # plt.plot(angles, peak(np.ones(np.shape(angles))*zPlot, angles, np.ones(np.shape(angles))*xPPlot), label='peak factor only')
            # #     plt.plot(zPlot, self.interpolator(zPlot, np.ones(np.shape(zPlot))*angle, np.ones(np.shape(zPlot))*xPPlot), label=str(xPPlot))
            # # plt.legend()
            # plt.show()

            # import matplotlib as mpl
            # mpl.use('TkAgg')
            # import matplotlib.pyplot as plt
            # plt.figure()
            # zValues = np.arange(-0.02, 0.02, 0.000001)
            # nonsmoothedTriangle = triangle(zValues)
            # # smoothedTriangle = triangle(zValues)+smoothCorner(zValues, -cuffWidth)+smoothCorner(zValues, cuffWidth);
            # # plt.plot(smoothedTriangle)
            # # plt.plot(np.diff(smoothedTriangle))
            # plt.plot(np.diff(np.diff(nonsmoothedTriangle)), label='smoothed')
            # # plt.plot(np.diff(np.diff(nonsmoothedTriangle)), label='not smoothed')
            # # plt.plot(np.diff(np.diff(smooth(nonsmoothedTriangle, 100))), label='numerically smoothed')
            # plt.legend()
            # plt.show()

        else:
            self.interpolator = interpolator

    def calculate_extracellular_potential(self, sourcePositions, sourceCurrents, receiverPositions):
        """

        :param sourcePositions: positions of current point sources
        :param sourceCurrents: currents of these point sources
        :param receiverPositions: positions the voltage is calculated for

        """

        def _interpolateStaticPotential(points):
            """
            Args:
                points: input from _compute_relative_positions_and_interpolate (this is the interpolation part)

            Returns: static potential that needs to be scaled with current

            """

            if self.method == 'z':
                zValues = points[2, :]
                # angle = 0  # arctan2(points[0,:], points[1,:])
                # xP = 0  # points[-1,:]
                return self.interpolator(zValues)  # triangularVoltage # smoothedVoltageStatic #

            elif self.method == 'z,xP,angle':
                zValues = points[2, :]
                angle = np.arctan2(points[0, :], points[1, :])
                xP = points[-1, :]
                return self.interpolator(zValues, angle, xP)  # triangularVoltage # smoothedVoltageStatic #

            else:
                raise KeyError




        # calculate LFP from membrane currents
        return compute_relative_positions_and_interpolate_fn_input(sourcePositions, sourceCurrents, receiverPositions,
                                                                   self.bundleGuide, _interpolateStaticPotential)

class precomputedFEM(ExtracellularPotentialMechanism):

    def __init__(self, bundleGuide, fieldName='noCuff1'):
        """Stores a precomputed voltage field and calculates the extracelluar potential caused by point sources. Used by :class:`recodingMechanismClass.RecordingMechanism`.

        :param bundleGuide: 3D nerve trajectory
        :param fieldName: string containing the name of the field to import; the location of field files needs to be specified somewhere (TODO)
        """

        # print '\nWhen using a recording mechanism based on a precomputed FEM model, make sure the electrodes are on a ' \
        #       'long (~1cm) straight part of the bundle guide.'

        # fieldDictArray = np.load(os.path.join('/Volumes/SANDISK/ComsolData/usedFields', fieldName, 'fieldDict.npy'))
        fieldDictArray = np.load(os.path.join('Fields', fieldName, 'fieldDict.npy'))
        self.FEMFieldDict = fieldDictArray[()]

        self.bundleGuide = bundleGuide

    def calculate_extracellular_potential(self, sourcePositions, sourceCurrents, receiverPositions):
        """

        :param sourcePositions: positions of current point sources
        :param sourceCurrents: currents of these point sources
        :param receiverPositions: positions the voltage is calculated for

        """

        # now interpolate from fieldImage
        positionToVoltageFn = lambda interpolationPoints: interpolateFromImage(self.FEMFieldDict, interpolationPoints, order=1)

        # calculate LFP from membrane currents
        return compute_relative_positions_and_interpolate_fn_input(sourcePositions, sourceCurrents, receiverPositions,
                                                            self.bundleGuide, positionToVoltageFn)

class homogeneous(ExtracellularPotentialMechanism):

    def __init__(self, sigma=1.):
        """The extracellular potential is calculated assuming a homogeneous outer medium with a constant, isotropic conductivity ``sigma``. Used by :class:`recodingMechanismClass.RecordingMechanism`.

        :param sigma: conductivity of the homogeneous tissue in S/m
        """

        self.sigma = sigma

    def calculate_extracellular_potential(self, sourcePositions, sourceCurrents, receiverPositions):
        """

        :param sourcePositions: positions of current point sources
        :param sourceCurrents: currents of these point sources
        :param receiverPositions: positions the voltage is calculated for

        """

        def _i_to_v_homogeneous(sourcePositions, sourceCurrents, receiverPositions, sigma=1., currentUnitSource=-9):
            """
            Idea and some implementation details from LFPy package

            Args:
                sourcePositions:
                sourceCurrents:
                receiverPositions:
                sigma:
                currentUnitSource:

            Returns:

            """

            # import matplotlib.pyplot as plt
            # plt.plot(sourceCurrents)
            # plt.show()

            nSourcePoints = np.shape(sourcePositions)[0]
            nReceiverPoints = np.shape(receiverPositions)[0]

            nTimePoints = len(sourceCurrents[:, 0])

            receiverPotentials = []
            for rInd in range(nReceiverPoints):
                receiverPosition = receiverPositions[rInd, :]

                r2 = (sourcePositions[:, 0] - receiverPosition[0]) ** 2 + (sourcePositions[:, 1] - receiverPosition[
                    1]) ** 2 + (sourcePositions[:, 2] - receiverPosition[2]) ** 2
                r = np.sqrt(r2)

                receiverPotential = 1 / (4 * np.pi * sigma) * np.dot(sourceCurrents.T, 1 / r)

                receiverPotentials.append(receiverPotential)

            receiverPotentials = np.array(receiverPotentials)

            return receiverPotentials

        return _i_to_v_homogeneous(sourcePositions, sourceCurrents, receiverPositions, sigma=self.sigma)





