from abc import ABCMeta, abstractmethod
import numpy as np
from createGeometry import random_perpendicular_vectors, rotation_matrix, length_from_coords
from nameSetters import get_directory_name
import os
import shutil
import silencer
import LFPy
import time

class RecordingMechanism(object):
    # __metaclass__ = ABCMeta

    def __init__(self, numberOfPoints, numberOfPoles, numberOfElectrodes, method='pointsource', sigma=0.3):

        if not numberOfPoles in [1,2]:
            raise Exception('Wrong number of poles, only 1 or 2 accepted. Received ' + str(numberOfPoles))

        self.numberOfPoints = numberOfPoints
        self.numberOfPoles = numberOfPoles
        self.numberOfElectrodes = numberOfElectrodes

        self.sigma = sigma
        self.method = method
        self.electrodeParameters = []

        self.electrodeDistances = []
        self.CAP = 0
        self.CAP_axonwise = []
        self.savePath = "" # where are the recordings stored?


    @abstractmethod
    def setup_recording_elec(self, bundleGuide, bundleLength):
        pass

    def load_one_axon(self, axonIndex):

        directory = self.savePath
        filename = "electrode_"+str(axonIndex)+".dat"

        electrodeData = np.loadtxt(os.path.join(directory, filename), unpack=True)

        return electrodeData


    def compute_overall_CAP(self):

        arrayd_axonwise_CAPs = np.array(self.CAP_axonwise)
        self.CAP = np.sum(arrayd_axonwise_CAPs, axis=0)


    def compute_single_axon_CAP(self, axon):
        """
        1. Calculate extracellular potential (LFP) for every electrode point defined in self.electrodeParameters
        2. Sum over points of single electrodes
        Args:
            axon:

        Returns: none

        """

        # 1. calculate LFP with LFPy from membrane currents

        # get the locations of electrodes, method of LFPy calculation and specific resistance
        electrodeParameters = self.electrodeParameters

        # shut down the output, always errors at the end because membrane current too high
        with silencer.nostdout():
            electrodes = LFPy.recextelectrode.RecExtElectrode(axon, **electrodeParameters)

            # calculate LFP by LFPy from membrane current
            electrodes.calc_lfp()

        # 2. sum points of individual electrodes up

        # crucial to define how iteration over electrodes goes
        monopolar = self.numberOfPoles == 1

        # get the local field potential
        LFP = electrodes.LFP

        # how long was the recording?
        length_t = np.shape(LFP)[1]

        # for every axon we have numberOfElectrodes electrodes
        CAP_axonwise = np.zeros((self.numberOfElectrodes, length_t))

        # The contactpoints that constitute one electrode contact (e.g. a cuff electrode ring) have to be recovered,
        # summed up together per recording location along the axon
        for i in range(self.numberOfElectrodes):
            if monopolar:
                contactPointIndices = range(self.numberOfPoints * i, self.numberOfPoints * (1 + i))
                meanOverContactPoints = np.mean(LFP[contactPointIndices, :], 0)
            else:
                contactPointIndicesPole1 = range(self.numberOfPoints * 2 * i, self.numberOfPoints * (1 + 2 * i))
                contactPointIndicesPole2 = range(self.numberOfPoints * (2 * i + 1),
                                                 self.numberOfPoints * (2 * (i + 1)))
                meanOverContactPoints = np.mean(
                    LFP[contactPointIndicesPole1, :] - LFP[contactPointIndicesPole2, :], 0)
                # sumOverContactPoints = np.sum(electrodeData[contactPointIndicesPole1, :] - electrodeData[contactPointIndicesPole2, :], 0)

            CAP_axonwise[i, :] = meanOverContactPoints

            if i == self.numberOfElectrodes - 1:
                self.CAP_axonwise.append(CAP_axonwise)




class RecCuff2D(RecordingMechanism):

    def __init__(self, radius, numberOfElectrodes = 1, numberOfPoints = 8, numberOfPoles=1, poleDistance=1000, positionMax=1, sigma=0.3, method='pointsource'):

        self.radius = radius # radius of the electrode ring (um)
        self.poleDistance = poleDistance # distance between two poles of the same bipolar electrode (um)
        self.positionMax = positionMax

        super(RecCuff2D,self).__init__(numberOfPoints, numberOfPoles, numberOfElectrodes, method, sigma)

    def setup_recording_elec(self, bundleGuide, bundleLength):

        # first find the bundle guide segment index that corresponds to the intendet bundle length (overlap for
        # myelinated axons gives longer bundle than specified by user
        bundleLengthIndex = np.shape(bundleGuide)[0]-1
        bundleLengthTemp = length_from_coords(bundleGuide)
        while bundleLengthTemp > bundleLength:
            bundleLengthIndex -= 1
            bundleLengthTemp = length_from_coords(bundleGuide[:bundleLengthIndex])

        # lastRecordedSegmentIndex = (np.shape(bundleGuide)[0]-1)*self.positionMax
        lastRecordedSegmentIndex = bundleLengthIndex*self.positionMax

        # distribute electrodes along bundle
        if self.numberOfElectrodes > 1:
            segmentIndices = np.linspace(lastRecordedSegmentIndex/self.numberOfElectrodes, lastRecordedSegmentIndex, self.numberOfElectrodes)
        else:
            segmentIndices = [lastRecordedSegmentIndex]

        # calculate absolute distances of electrode along the bundle from orign of bundle
        for i in range(self.numberOfElectrodes):
            self.electrodeDistances.append(np.floor(length_from_coords(bundleGuide[:segmentIndices[i]])))

        # variable to save points of electrode
        electrodePositions = np.array([]).reshape(0,3)

        for i in range(self.numberOfElectrodes):
            segmentNumber = int(segmentIndices[i])
            # segmentNumber = floor(np.shape(bundleGuide)[0]/numberOfElectrodes)*(i+1) - 1

            segmentStartingPos = bundleGuide[segmentNumber - 1,:]
            segmentEndPos = bundleGuide[segmentNumber,:]

            segmentMiddle = (segmentStartingPos + segmentEndPos)/2
            segmentOrientation = segmentEndPos - segmentStartingPos
            segmentOrientation = segmentOrientation/np.linalg.norm(segmentOrientation)

            # get one random orthogonal vector
            orthogonalVector = random_perpendicular_vectors(segmentOrientation)[0,:]

            electrodePositionsOneElectrode = np.array([]).reshape(0,3)

            # loop to generate one ring
            for j in range(self.numberOfPoints):
                # generate the coordinates for one ring for the first pole of the electrode
                pointPosition = np.dot(rotation_matrix(segmentOrientation, 2*np.pi/self.numberOfPoints*j),(orthogonalVector*self.radius)) + segmentMiddle

                # append it to the list of coordinates for this pole
                electrodePositionsOneElectrode = np.vstack([electrodePositionsOneElectrode, pointPosition])

            # if the electrodes are bipolar
            if self.numberOfPoles == 2:
                electrodePositionsPole2 = electrodePositionsOneElectrode + np.tile(segmentOrientation*self.poleDistance, (np.shape(electrodePositionsOneElectrode)[0],1))
                electrodePositionsOneElectrode = np.vstack([electrodePositionsOneElectrode, electrodePositionsPole2])

            electrodePositions = np.vstack([electrodePositions, electrodePositionsOneElectrode])


        electrodeParameters = {
                    'sigma' : self.sigma,
                    'x' : electrodePositions[:,0],  #Coordinates of electrode contacts
                    'y' : electrodePositions[:,1],
                    'z' : electrodePositions[:,2],
                    # 'n' : 20,
                    # 'r' : 10,
                    # 'N' : N,
                    'method': self.method, #or "linesource"
                }

        self.electrodeParameters = electrodeParameters


class RecCuff3D(RecordingMechanism):

    def __init__(self, radius, numberOfElectrodes = 1, width = 20000, numberOfPoles=1, poleDistance=1000, pointsPerRing=8, distanceOfRings=400, positionMax=1, sigma=0.3, method='pointsource'):

        self.radius = radius # radius of the electrode ring (um)
        self.width = width # dimension of electrode along bundle
        self.poleDistance = poleDistance # distance between two poles of the same bipolar electrode (um)
        self.positionMax = positionMax

        # self.pointsPerMicroMeter = pointsPerMicroMeter
        # self.pointsPerPointRing =  max(1, np.floor(pointsPerMicroMeter*2*np.pi*radius))
        # self.numberOfPointRings = max(1, np.floor(pointsPerMicroMeter*width))
        self.pointsPerPointRing =  pointsPerRing
        self.distanceOfRings = distanceOfRings
        self.numberOfPointRings = max(1, np.floor(width/distanceOfRings))

        # numberOfPoints = np.floor(pointsPerMicroMeter*width)*np.floor(pointsPerMicroMeter*2*np.pi) # how many monopolar electrodes are used to approximate the electrode
        numberOfPoints = int(self.pointsPerPointRing*self.numberOfPointRings)

        super(RecCuff3D,self).__init__(numberOfPoints, numberOfPoles, numberOfElectrodes, method, sigma)

    def setup_recording_elec(self, bundleGuide, bundleLength):

        # first find the bundle guide segment index that corresponds to the intendet bundle length (overlap for
        # myelinated axons gives longer bundle than specified by user
        bundleLengthIndex = np.shape(bundleGuide)[0]-1
        bundleLengthTemp = length_from_coords(bundleGuide)
        while bundleLengthTemp > bundleLength:
            bundleLengthIndex -= 1
            bundleLengthTemp = length_from_coords(bundleGuide[:bundleLengthIndex])

        # lastRecordedSegmentIndex = (np.shape(bundleGuide)[0]-1)*self.positionMax
        lastRecordedSegmentIndex = bundleLengthIndex*self.positionMax

        if self.numberOfElectrodes > 1:
            segmentIndices = np.linspace(lastRecordedSegmentIndex/self.numberOfElectrodes, lastRecordedSegmentIndex, self.numberOfElectrodes)
        else:
            segmentIndices = [lastRecordedSegmentIndex]

        for i in range(self.numberOfElectrodes):
            self.electrodeDistances.append(np.floor(length_from_coords(bundleGuide[:segmentIndices[i]])))

        electrodePositions = np.array([]).reshape(0,3)

        for i in range(self.numberOfElectrodes):
            segmentNumber = int(segmentIndices[i])
            # segmentNumber = floor(np.shape(bundleGuide)[0]/numberOfElectrodes)*(i+1) - 1

            segmentStartingPos = bundleGuide[segmentNumber - 1,:]
            segmentEndPos = bundleGuide[segmentNumber,:]

            segmentMiddle = (segmentStartingPos + segmentEndPos)/2
            segmentOrientation = segmentEndPos - segmentStartingPos
            segmentOrientation = segmentOrientation/np.linalg.norm(segmentOrientation)

            # get one random orthogonal vector
            orthogonalVector = random_perpendicular_vectors(segmentOrientation)[0,:]

            electrodePositionsOneElectrode = np.array([]).reshape(0,3)

            electrodePositionsOneRing = np.array([]).reshape(0,3) # one ring, gets shifted by density to create surface
            # loop to generate one ring
            for j in range(int(self.pointsPerPointRing)):
                # generate the coordinates for one ring for the first pole of the electrode
                pointPosition = np.dot(rotation_matrix(segmentOrientation, 2*np.pi/self.pointsPerPointRing*j),(orthogonalVector*self.radius)) + segmentMiddle

                # append it to the list of coordinates for this pole
                electrodePositionsOneRing = np.vstack([electrodePositionsOneRing, pointPosition])

            # shift the created ring long bundle and add to electrode coordinates to create the surface
            for k in range(int(self.numberOfPointRings)):

                shiftedRing = electrodePositionsOneRing + np.tile(segmentOrientation*self.distanceOfRings*k,
                                                                  (np.shape(electrodePositionsOneRing)[0],1))
                # append shifted ring points to the list of coordinates for this pole
                electrodePositionsOneElectrode = np.vstack([electrodePositionsOneElectrode, shiftedRing])

            # if the electrodes are bipolar
            if self.numberOfPoles == 2:
                electrodePositionsPole2 = electrodePositionsOneElectrode + np.tile(segmentOrientation*self.poleDistance, (np.shape(electrodePositionsOneElectrode)[0],1))
                electrodePositionsOneElectrode = np.vstack([electrodePositionsOneElectrode, electrodePositionsPole2])

            electrodePositions = np.vstack([electrodePositions, electrodePositionsOneElectrode])


        electrodeParameters = {
                    'sigma' : self.sigma,
                    'x' : electrodePositions[:,0],  #Coordinates of electrode contacts
                    'y' : electrodePositions[:,1],
                    'z' : electrodePositions[:,2],
                    # 'n' : 20,
                    # 'r' : 10,
                    # 'N' : N,
                    'method': self.method, #or "linesource"
                }

        self.electrodeParameters = electrodeParameters


# class CuffElectrode3D(RecordingMechanism):
#     def __init__(self, radius, width = 15, poles=1, poleDistance=20, pointDensity=0.1):
#
#         self.radius = radius # radius of the electrode ring (um)
#         self.width = width # width in axial direction along the bundle (um)
#         self.poles = poles # mono- or bipolar
#         self.poleDistance = poleDistance # distance between two poles of the same bipolar electrode (um)
#         self.pointDensity = pointDensity # density of points that approximate a surface (1/um)
#
#         numberOfElectrodes = radius*width*2*np.pi*pointDensity
#
#         super(CuffElectrode3D,self).__init__(numberOfElectrodes)

# elec = CuffElectrodes()
#
# print elec