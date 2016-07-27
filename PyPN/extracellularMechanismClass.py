from abc import abstractmethod, ABCMeta
import numpy as np
import os
import silencer
import LFPy
import time
from extracellularBackend import *

class ExtracellularPotentialMechanism(object):
    __metaclass__ = ABCMeta

    @abstractmethod

    def calculate_extracellular_potential(self, sourcePositions, sourceCurrents, receiverPositions):
        pass

class precomputedFEM(ExtracellularPotentialMechanism):

    def __init__(self, bundleGuide, fieldName='noCuff1'):

        # todo: input params that allow more degrees of freedom for the FEM model

        print '\nWhen using a recording mechanism based on a precomputed FEM model, make sure the electrodes are on a ' \
              'long (~1cm) straight part of the bundle guide.'

        fieldDictArray = np.load(
            os.path.join('/media/carl/4ECC-1C44/ComsolData/usedFields', fieldName, 'fieldDict.npy'))
        self.FEMFieldDict = fieldDictArray[()]

        self.bundleGuide = bundleGuide

    # def calculate_LFP(self, axon, electrodePositions): # sourcePositions, sourceCurrents,
    def calculate_extracellular_potential(self, sourcePositions, sourceCurrents, receiverPositions):

        # calculate LFP from membrane currents
        return compute_relative_positions_and_interpolate(sourcePositions, sourceCurrents, receiverPositions, self.FEMFieldDict, self.bundleGuide)



class homogeneous(ExtracellularPotentialMechanism):

    def __init__(self, sigma=1.):

        self.sigma = sigma

    def calculate_extracellular_potential(self, sourcePositions, sourceCurrents, receiverPositions):

        # np.vstack([axon.xmid, axon.ymid, axon.zmid]).T, axon.imem, electrodePositions
        return i_to_v_homogeneous(sourcePositions, sourceCurrents, receiverPositions, sigma=self.sigma)

    # def calculate_LFP(self, axon, electrodePositions):
    #     """
    #     Calculate extracellular potential (LFP) for every electrode point defined in self.electrodeParameters
    #     Args:
    #         axon:
    #
    #     Returns: LFP
    #
    #     """
    #
    #     #calculate LFP with LFPy from membrane currents
    #
    #     # put the locations of electrodes, the method and sigma in the right form for LFPy
    #     electrodeParameters = {
    #         'sigma': self.sigma,
    #         'x': electrodePositions[:, 0],  # Coordinates of electrode contacts
    #         'y': electrodePositions[:, 1],
    #         'z': electrodePositions[:, 2],
    #         'method': self.method,  # 'pointsource' or "linesource"
    #     }
    #
    #     # shut down the output, always errors at the end because membrane current too high
    #     with silencer.nostdout():
    #         electrodes = LFPy.recextelectrode.RecExtElectrode(axon, **electrodeParameters)
    #
    #         # calculate LFP by LFPy from membrane current
    #         electrodes.calc_lfp()
    #
    #
    #     # get the local field potential
    #     return electrodes.LFP






