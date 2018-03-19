from abc import abstractmethod, ABCMeta
import numpy as np
import os
import silencer
import time
from extracellularBackend import *
from scipy.interpolate import interp1d
import analyticFnGen

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

            self.interpolator = analyticFnGen.idealizedCuff(cuffWidthIn=0.02)

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
                angle = np.arctan2(points[1, :], points[0, :])
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





