# from abc import ABCMeta, abstractmethod
import numpy as np
from createGeometry import random_perpendicular_vectors, rotation_matrix, length_from_coords
from nameSetters import get_directory_name
import os
import shutil

class RecordingMechanism(object):
    # __metaclass__ = ABCMeta

    def __init__(self, numberOfPoints, numberOfPoles, numberOfElectrodes, method='pointsource', sigma=0.3):

        self.numberOfPoints = numberOfPoints
        self.numberOfPoles = numberOfPoles
        self.numberOfElectrodes = numberOfElectrodes

        self.sigma = sigma
        self.method = method
        self.electrodeParameters = []

        self.electrodeDistances = []
        self.CAP = 0
        self.CAP_axonwise = 0
        self.savePath = "" # where are the recordings stored?



    def setup_recording_elec(self, bundleGuide, bundleLength):
        pass

    def load_one_axon(self, axonIndex):

        directory = self.savePath
        filename = "electrode_"+str(axonIndex)+".dat"

        electrodeData = np.loadtxt(os.path.join(directory, filename), unpack=True)

        return electrodeData

    def compute_CAP_from_files(self):

        monopolar = self.numberOfPoles == 1

        # how long was the recording?
        electrodeData = self.load_one_axon(0)
        length_t = np.shape(electrodeData)[1]

        # variable to save the sum over all axons
        self.CAP = np.zeros((self.numberOfElectrodes, length_t))

        # variable to save the extracellular signal from each cell separately, at the last electrode position.
        self.CAP_axonwise = np.array([]).reshape(0,length_t) # np.zeros((self.numberOfAxons, length_t))

        # load the recordings for every axon one by one and add them.
        axonIndex = 0
        while True:
            try:
                electrodeData = self.load_one_axon(axonIndex)
            except:
                break

            # The contactpoints that constitute one cuff electrode ring have to be recovered, summed up together per
            # recording location along the axon
            for i in range(self.numberOfElectrodes):
                if monopolar:
                    contactPointIndices = range(self.numberOfPoints*i, self.numberOfPoints*(1+i))
                    sumOverContactPoints = np.mean(electrodeData[contactPointIndices, :], 0)
                    # sumOverContactPoints = np.sum(electrodeData[contactPointIndices, :], 0)
                else:
                    contactPointIndicesPole1 = range(self.numberOfPoints*2*i, self.numberOfPoints*(1+2*i))
                    contactPointIndicesPole2 = range(self.numberOfPoints*(2*i+1), self.numberOfPoints*(2*(i+1)))
                    sumOverContactPoints = np.mean(electrodeData[contactPointIndicesPole1, :] - electrodeData[contactPointIndicesPole2, :], 0)
                    # sumOverContactPoints = np.sum(electrodeData[contactPointIndicesPole1, :] - electrodeData[contactPointIndicesPole2, :], 0)

                self.CAP[i,:] = self.CAP[i,:] +  sumOverContactPoints

                if i == self.numberOfElectrodes-1:
                    # self.AP_axonwise[axonIndex,:] = sumOverContactPoints
                    self.CAP_axonwise = np.vstack([self.CAP_axonwise, sumOverContactPoints])

            axonIndex += 1




class CuffElectrode2D(RecordingMechanism):

    def __init__(self, radius, numberOfElectrodes = 1, numberOfPoints = 8, numberOfPoles=1, poleDistance=20, positionMax=1, sigma=0.3, method='pointsource'):

        self.radius = radius # radius of the electrode ring (um)
        self.numberOfPoles = numberOfPoles # mono- or bipolar

        if not numberOfPoles in [1,2]:
            raise Exception('Wrong number of poles, only 1 or 2 accepted. Received ' + str(numberOfPoles))

        self.poleDistance = poleDistance # distance between two poles of the same bipolar electrode (um)
        self.numberOfPoints = numberOfPoints # how many monopolar electrodes are used to approximate the electrode
        self.positionMax = positionMax

        super(CuffElectrode2D,self).__init__(numberOfPoints, numberOfPoles, numberOfElectrodes, method)

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

    def calculate_CAP(self):
        print 'bliblub'


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