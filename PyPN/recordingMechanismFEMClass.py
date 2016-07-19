from abc import abstractmethod
import numpy as np
from scipy import ndimage
from createGeometry import random_perpendicular_vectors, rotation_matrix, length_from_coords
from nameSetters import get_directory_name
import os
import takeTime
import shutil
import silencer
import LFPy
import time
import cPickle as pickle

import matplotlib.pyplot as plt
#TODO: delete, only for testing
import matplotlib.cm as cm
import matplotlib.colors as colors

class RecordingMechanismFEM(object):

    def __init__(self, numberOfPoints, numberOfPoles, numberOfElectrodes):

        if not numberOfPoles in [1,2]:
            raise Exception('Wrong number of poles, only 1 or 2 accepted. Received ' + str(numberOfPoles))

        print '\nWhen using a recording mechanism based on a precomputed FEM model, make sure the electrodes are on a ' \
              'long (~1cm) straight part of the bundle guide.'

        self.numberOfPoints = numberOfPoints
        self.numberOfPoles = numberOfPoles
        self.numberOfElectrodes = numberOfElectrodes
        self.electrodePositions = np.array([])

        self.bundleGuide = np.array([])

        self.FEMFieldLocation = []
        # TODO: input folder and axonXs
        # folder =  '/media/carl/UUI/Comsol/data/noCuffFinerGrid/z0.003/Gauss_Points_Order3' # '/media/carl/4ECC-1C44/ComsolData/thickerEndoneurium' # '/media/carl/4ECC-1C44/ComsolData/cuff_0.5mm' #
        # folder = '/media/carl/4ECC-1C44/ComsolData/noCuffFiner/z0.003_100/Lagrange_Smoothing_Finer' # '/media/carl/4ECC-1C44/ComsolData/thickerEndoneurium' # '/media/carl/4ECC-1C44/ComsolData/noCuffFiner/z0.001_500' # '/media/carl/4ECC-1C44/ComsolData/noCuffFiner/z0.001_500' #
        folder  = '/media/carl/4ECC-1C44/ComsolData/noCuffFiner/z0.03_1000,x0.0015_100,y_asym'
        # folder = '/media/carl/4ECC-1C44/ComsolData/cuffFiner'

        axonXs = [0, 180]
        # self.FEMFieldDict = pickle.load(open(os.path.join(folder, 'numpy', 'fieldDict.pickle'), 'rb')) # self.load_field(folder, axonXs)
        fieldDictArray = np.load(os.path.join(folder, 'numpy', 'fieldDict.npy'))
        self.FEMFieldDict = fieldDictArray[()]

        self.electrodeDistances = []
        self.CAP = 0
        self.CAP_axonwise = []
        self.savePath = "" # where are the recordings stored?



    def load_field(self, folder, axonXs):

        # get file names
        filenames = [f for f in sorted(os.listdir(folder)) if os.path.isfile(os.path.join(folder, f))]

        axonXSteps = len(axonXs)
        assert axonXSteps == len(filenames)

        # load each field (different axon positions)
        fields = []
        for filename in filenames:
            fields.append(np.loadtxt(os.path.join(folder, filename)))

        print 'loaded field'

        # get coordinates (should be equal for all field files, otherwise nothing works)
        x = fields[0][:, 0]
        y = fields[0][:, 1]
        z = fields[0][:, 2]

        # sort by coordinate values, x changing fastest, z slowest
        orderIndices = np.lexsort((x, y, z))
        x = x[orderIndices]
        y = y[orderIndices]
        z = z[orderIndices]

        # get coordinate values
        xValues = np.unique(x)
        yValues = np.unique(y)
        zValues = np.unique(z)

        # get number of steps
        xSteps = len(xValues)
        ySteps = len(yValues)
        zSteps = len(zValues)

        # voltages are different for each field
        voltages = []
        for i in range(axonXSteps):
            v = fields[i][:, 3]
            v = v[orderIndices]  # order voltages as well
            voltages.append(v)

        # transform data to 3D-field with integer indices replacing actual coordinate values
        fieldImage = np.zeros([xSteps, ySteps, zSteps, axonXSteps])

        for axonXInd in range(axonXSteps):
            for xInd in range(xSteps):
                for yInd in range(ySteps):
                    for zInd in range(zSteps):
                        vIndexCalc = xInd + xSteps * (yInd + zInd * ySteps)
                        fieldImage[xInd, yInd, zInd, axonXInd] = voltages[axonXInd][vIndexCalc]

        fieldDict = {'fieldImage': fieldImage,
                     'x': xValues,
                     'y': yValues,
                     'z': zValues,
                     'axonX': axonXs}

        return fieldDict

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
        1. Calculate extracellular potential for every electrode point defined in self.electrodeParameters
        2. Sum over points of single electrodes
        Args:
            axon:

        Returns: none

        """

        def _getImageCoords(fieldDict, points):

            xValues = fieldDict['x']
            yValues = fieldDict['y']
            zValues = fieldDict['z']
            axonXValues = fieldDict['axonX']

            # assume equidistant original points

            xMin = min(xValues)
            xMax = max(xValues)
            xNum = len(xValues)

            yMin = min(yValues)
            yMax = max(yValues)
            yNum = len(yValues)

            zMin = min(zValues)
            zMax = max(zValues)
            zNum = len(zValues)

            axonXMin = min(axonXValues)
            axonXMax = max(axonXValues)
            axonXNum = len(axonXValues)

            points = np.array(points)

            if len(points.shape) > 1:
                if points.shape[1] > 4:
                    points = np.transpose(points)
                xCoords = np.add(points[:, 0], -xMin) / (xMax - xMin) * (xNum - 1)
                yCoords = np.add(points[:, 1], -yMin) / (yMax - yMin) * (yNum - 1)
                zCoords = np.add(points[:, 2], -zMin) / (zMax - zMin) * (zNum - 1)
                xAxonCoords = np.add(points[:, 3], -axonXMin) / (axonXMax - axonXMin) * (axonXNum - 1)
            else:
                xCoords = (points[0] - xMin) / (xMax - xMin) * (xNum - 1)
                yCoords = (points[1] - yMin) / (yMax - yMin) * (yNum - 1)
                zCoords = (points[2] - zMin) / (zMax - zMin) * (zNum - 1)
                xAxonCoords = (points[3] - axonXMin) / (axonXMax - axonXMin) * (axonXNum - 1)

            zCoords = np.abs(
                zCoords) # in the input FEM field, we only take one side of the z-value range thanks to symmetry
            yCoords = np.abs(yCoords) # same for y-Coordinates

            return np.vstack([xCoords, yCoords, zCoords, xAxonCoords])

        def _interpolateFromImage(fieldDict, points, order=3):

            # first transform coordinates in points into position coordinates
            imageCoords = _getImageCoords(fieldDict, points)

            # then with new coords to the interpolation
            return ndimage.map_coordinates(fieldDict['fieldImage'], imageCoords, order=order)

        def _compute_relative_positions_and_interpolate(axon):

            """

            Structure as follows

            for bs in bundleSegments:
                for a in axonsOfBundleSegment:
                    calculate axonDir2DNorms and axonPositionParallelToBundle
                for e in electrodes calculate:
                    for a in axonsOfBundleSegment:
                          calculate electrode position and interpolate potential

            Args:
                axon: needed for segment positions and membrane currents

            Returns:

            """

            axonSegCoords = np.vstack([axon.xmid, axon.ymid, axon.zmid])
            nseg = len(axon.xmid)

            bundleCoords = self.bundleGuide[:,0:3]

            # first go through all bundle segments, find the associated axon segments and calculate the needed
            # quantities for all axon segments
            elecPotentials = np.zeros((self.electrodePositions.shape[0], np.shape(axon.imem)[1]))

            # # TODO: delete this again, only for testing
            # lastElectrodeSignal = []
            # f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
            # jet = plt.get_cmap('jet')
            # cNorm = colors.Normalize(vmin=0, vmax=int(nseg) - 1)
            # scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

            t0 = time.time()

            bundleSegInd = 1
            axonSegInd = 0
            axonFinishedFlag = False
            while bundleSegInd < bundleCoords.shape[0] and not axonFinishedFlag: # bundleSegInd points to endpoint of current bundle segment

                bundleSegEndPoint = bundleCoords[bundleSegInd, :]
                bundleSegStartPoint = bundleCoords[bundleSegInd - 1, :]

                # find the normal to the terminating surface of the bundle guide segment
                dir0 = bundleSegEndPoint - bundleSegStartPoint
                dir0norm = dir0/np.linalg.norm(dir0)
                if bundleSegInd < bundleCoords.shape[0]-1:
                    dir1 = bundleCoords[bundleSegInd + 1, :] - bundleSegEndPoint
                    n = (dir0 + dir1)/2
                else:
                    n = dir0

                # # TODO: maybe this introduces problems?
                # randomAxonXDir = random_perpendicular_vectors(dir0norm)[0,:]

                # process axons associated with the current bundle segment
                # and store axon segment properties in these lists
                axonDir2Ds = []
                axonXDists = []
                axonZs = []

                axonSegPos = axonSegCoords[:, axonSegInd] # why oh why, no do-while
                while np.inner(axonSegPos - bundleSegEndPoint, n) < 0: # while in bundle guide segment

                    axonSegPos = axonSegCoords[:, axonSegInd]
                    axonPosPar = np.inner(axonSegPos - bundleSegStartPoint,
                                          dir0norm)  # compontent parallel to bundle direction

                    # normal from bundle guide to axon position
                    axonDir2D = axonSegPos - (bundleSegStartPoint + axonPosPar * dir0norm)
                    axonXDist = np.linalg.norm(axonDir2D)
                    if not axonXDist == 0:
                        axonDir2DNorm = axonDir2D / axonXDist
                    else:
                        axonDir2DNorm = -1

                    # save all computed values of axon into these lists, they are used when iterating over electrodes
                    axonXDists.append(axonXDist)
                    axonDir2Ds.append(axonDir2DNorm)
                    axonZs.append(axonPosPar)

                    axonSegInd += 1

                    if axonSegInd >= nseg:
                        axonFinishedFlag=True
                        break

                # now process electrodes

                # compontent parallel to bundle direction
                bundleSegStartPointTiled = np.tile(bundleSegStartPoint, (self.electrodePositions.shape[0], 1))
                elecPosPar = np.inner(self.electrodePositions - bundleSegStartPointTiled, dir0norm)

                # electrode shifted to origin
                electrodeVector = self.electrodePositions - (bundleSegStartPointTiled + np.tile(dir0norm,
                                    (self.electrodePositions.shape[0], 1)) * elecPosPar[:, np.newaxis])

                # TODO: take only axons that are inside FEM field size?
                # TODO: (check how long it is, maybe vary the length)
                for axonSegIndLoc in range(len(axonDir2Ds)):
                    axonDir2D = axonDir2Ds[axonSegIndLoc]

                    if isinstance(axonDir2D, int): # if the axon segment lies on the bundle middle exactly
                        elecX = np.ones(self.electrodePositions.shape[0])*np.linalg.norm(electrodeVector[0,:])
                        elecY = np.zeros(self.electrodePositions.shape[0])
                    else:
                        axonDir2DTiled = np.tile(axonDir2D, (electrodeVector.shape[0], 1))

                        # electrode coordinates projected onto new base vectors
                        elecX = np.inner(electrodeVector, axonDir2D)

                        elecYVec = electrodeVector - axonDir2DTiled * elecX[:, np.newaxis]
                        elecY = np.linalg.norm(elecYVec, axis=1)

                    axonZ = axonZs[axonSegIndLoc]
                    elecZ = np.add(elecPosPar, -axonZ)

                    axonXDist = np.tile(axonXDists[axonSegIndLoc],(1, len(elecZ)))

                    # # TODO: remove this again, only testing
                    # elecX = np.zeros(self.electrodePositions.shape[0])
                    # elecY = np.zeros(self.electrodePositions.shape[0])

                    interpolationPoints = np.vstack([elecX, elecY, elecZ, axonXDist])
                    interpolationPoints = np.divide(interpolationPoints,1000000)  # from um to m TODO: numercal problem?

                    # now interpolate from fieldImage
                    elecPotTempStatic = _interpolateFromImage(self.FEMFieldDict, interpolationPoints, order=1)

                    # # TODO: replace this by fitted function, just for testing here
                    # def func1(x, a, b, c, d):
                    #     return a * (x + b) ** (-c) + d
                    # elecPotTempStatic = func1(elecZ, 10 ** (-2), 0.0005, 2, 0)

                    imemAxonSegInd = axonSegInd-len(axonDir2Ds)+axonSegIndLoc

                    elecPotTemp = np.outer(elecPotTempStatic, axon.imem[imemAxonSegInd,:])*1000 # COMSOL gave V, we need mV

                    # # if np.mod(imemAxonSegInd,50) == 0:
                    # #     # plt.plot(axon.imem[imemAxonSegInd,50:], label='Current') # /np.max(axon.imem[imemAxonSegInd,:])
                    # #     plt.plot(np.transpose(elecPotTemp[-1,50:]), label='Pot') # /np.max(elecPotTemp[-1, :])
                    # ax1.plot(np.transpose(elecPotTemp[-1, 0:]), color=scalarMap.to_rgba(imemAxonSegInd))  # /np.max(elecPotTemp[-1, :])
                    #
                    # lastElectrodeSignal.append(np.squeeze(elecPotTemp[-1, 0:]))

                    # add contributions
                    elecPotentials = np.add(elecPotentials, elecPotTemp)

                bundleSegInd += 1

            # lfp = np.zeros(len(np.squeeze(axon.imem[0,:])))
            # lfps = []
            # for axonSegInd in range(nseg):
            #     axonSegPos = axonSegCoords[:, axonSegInd]
            #     lastElecPos = self.electrodePositions[-1,:]
            #
            #     r = np.linalg.norm(axonSegPos - lastElecPos)
            #
            #     lfpTemp = 1/(4*np.pi*2*r)*axon.imem[axonSegInd,:]
            #
            #     lfps.append(lfpTemp)
            #
            #     lfp = lfp + lfpTemp



            # ax2.plot(np.transpose(np.array(lfps)))
            #
            #
            #
            # lastElectrodeSignalArray = np.array(lastElectrodeSignal)
            # np.save('/media/carl/4ECC-1C44/PyPN/FEM_CAPs/z0.003m_100steps_longer.npy', lastElectrodeSignalArray)
            # # ax3.plot(np.sum(lastElectrodeSignalArray[0:10,:], axis=0), linewidth=2, label='until 10')
            # # ax3.plot(np.sum(lastElectrodeSignalArray[0:30, :], axis=0), linewidth=2, label='until 30')
            # # ax3.plot(np.sum(lastElectrodeSignalArray[0:50, :], axis=0), linewidth=2, label='until 50')
            # # ax3.plot(np.sum(lastElectrodeSignalArray[0:80, :], axis=0), linewidth=2, label='until 80')
            # # ax3.plot(np.sum(lastElectrodeSignalArray[0:120, :], axis=0), linewidth=2, label='until 120')
            # # ax3.plot(np.sum(lastElectrodeSignalArray[0:200, :], axis=0), linewidth=2, label='until 200')
            # cNorm = colors.Normalize(vmin=0, vmax=10)
            # scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)
            # ind = 0
            # step = 10
            # plotInd = 0
            # while ind + step - 1 < lastElectrodeSignalArray.shape[0]:
            #     # ax3.plot(np.sum(lastElectrodeSignalArray[ind:(ind+step), :], axis=0), linewidth=2, label='['+str(ind)+','+str(ind+step-1)+']', color=scalarMap.to_rgba(plotInd))
            #     ax3.plot(np.sum(lastElectrodeSignalArray[0:(ind + step), :], axis=0), linewidth=2,
            #              label='[' + str(0) + ',' + str(ind + step - 1) + ']', color=scalarMap.to_rgba(plotInd))
            #     ind += step
            #     plotInd += 1
            # ax3.plot(np.sum(lastElectrodeSignalArray, axis=0), linewidth=2, label='all', color='red')
            #
            # ax3.plot(lfp, linewidth=2, label='analytical')
            # plt.legend()
            # plt.show()


            print 'finished all bundle segments in %i seconds.' % (time.time() - t0)

            return elecPotentials


        # 1. calculate LFP with LFPy from membrane currents

        LFP = _compute_relative_positions_and_interpolate(axon)

        # 2. sum points of individual electrodes up

        # crucial to define how iteration over electrodes goes
        monopolar = self.numberOfPoles == 1

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




class FEMRecCuff2D(RecordingMechanismFEM):

    def __init__(self, radius, numberOfElectrodes = 1, numberOfPoints = 8, numberOfPoles=1, poleDistance=1000, positionMax=1):

        self.radius = radius # radius of the electrode ring (um)
        self.poleDistance = poleDistance # distance between two poles of the same bipolar electrode (um)
        self.positionMax = positionMax

        super(FEMRecCuff2D, self).__init__(numberOfPoints, numberOfPoles, numberOfElectrodes)

    def setup_recording_elec(self, bundleGuide, bundleLength):

        bundleGuide = bundleGuide[:, 0:3]

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


        # electrodeParameters = {
        #             'sigma' : self.sigma,
        #             'x' : electrodePositions[:,0],  #Coordinates of electrode contacts
        #             'y' : electrodePositions[:,1],
        #             'z' : electrodePositions[:,2],
        #             # 'n' : 20,
        #             # 'r' : 10,
        #             # 'N' : N,
        #             'method': self.method, #or "linesource"
        #         }

        self.electrodePositions = electrodePositions


# class RecCuff3D(RecordingMechanism):
#
#     def __init__(self, radius, numberOfElectrodes = 1, width = 20000, numberOfPoles=1, poleDistance=1000, pointsPerRing=8, distanceOfRings=400, positionMax=1, sigma=0.3, method='pointsource'):
#
#         self.radius = radius # radius of the electrode ring (um)
#         self.width = width # dimension of electrode along bundle
#         self.poleDistance = poleDistance # distance between two poles of the same bipolar electrode (um)
#         self.positionMax = positionMax
#
#         # self.pointsPerMicroMeter = pointsPerMicroMeter
#         # self.pointsPerPointRing =  max(1, np.floor(pointsPerMicroMeter*2*np.pi*radius))
#         # self.numberOfPointRings = max(1, np.floor(pointsPerMicroMeter*width))
#         self.pointsPerPointRing =  pointsPerRing
#         self.distanceOfRings = distanceOfRings
#         self.numberOfPointRings = max(1, np.floor(width/distanceOfRings))
#
#         # numberOfPoints = np.floor(pointsPerMicroMeter*width)*np.floor(pointsPerMicroMeter*2*np.pi) # how many monopolar electrodes are used to approximate the electrode
#         numberOfPoints = int(self.pointsPerPointRing*self.numberOfPointRings)
#
#         super(RecCuff3D,self).__init__(numberOfPoints, numberOfPoles, numberOfElectrodes, method, sigma)
#
#     def setup_recording_elec(self, bundleGuide, bundleLength):
#
#         bundleGuide = bundleGuide[:, 0:3]
#
#         # first find the bundle guide segment index that corresponds to the intendet bundle length (overlap for
#         # myelinated axons gives longer bundle than specified by user
#         bundleLengthIndex = np.shape(bundleGuide)[0]-1
#         bundleLengthTemp = length_from_coords(bundleGuide)
#         while bundleLengthTemp > bundleLength:
#             bundleLengthIndex -= 1
#             bundleLengthTemp = length_from_coords(bundleGuide[:bundleLengthIndex])
#
#         # lastRecordedSegmentIndex = (np.shape(bundleGuide)[0]-1)*self.positionMax
#         lastRecordedSegmentIndex = bundleLengthIndex*self.positionMax
#
#         if self.numberOfElectrodes > 1:
#             segmentIndices = np.linspace(lastRecordedSegmentIndex/self.numberOfElectrodes, lastRecordedSegmentIndex, self.numberOfElectrodes)
#         else:
#             segmentIndices = [lastRecordedSegmentIndex]
#
#         for i in range(self.numberOfElectrodes):
#             self.electrodeDistances.append(np.floor(length_from_coords(bundleGuide[:segmentIndices[i]])))
#
#         electrodePositions = np.array([]).reshape(0,3)
#
#         for i in range(self.numberOfElectrodes):
#             segmentNumber = int(segmentIndices[i])
#             # segmentNumber = floor(np.shape(bundleGuide)[0]/numberOfElectrodes)*(i+1) - 1
#
#             segmentStartingPos = bundleGuide[segmentNumber - 1,:]
#             segmentEndPos = bundleGuide[segmentNumber,:]
#
#             segmentMiddle = (segmentStartingPos + segmentEndPos)/2
#             segmentOrientation = segmentEndPos - segmentStartingPos
#             segmentOrientation = segmentOrientation/np.linalg.norm(segmentOrientation)
#
#             # get one random orthogonal vector
#             orthogonalVector = random_perpendicular_vectors(segmentOrientation)[0,:]
#
#             electrodePositionsOneElectrode = np.array([]).reshape(0,3)
#
#             electrodePositionsOneRing = np.array([]).reshape(0,3) # one ring, gets shifted by density to create surface
#             # loop to generate one ring
#             for j in range(int(self.pointsPerPointRing)):
#                 # generate the coordinates for one ring for the first pole of the electrode
#                 pointPosition = np.dot(rotation_matrix(segmentOrientation, 2*np.pi/self.pointsPerPointRing*j),(orthogonalVector*self.radius)) + segmentMiddle
#
#                 # append it to the list of coordinates for this pole
#                 electrodePositionsOneRing = np.vstack([electrodePositionsOneRing, pointPosition])
#
#             # shift the created ring long bundle and add to electrode coordinates to create the surface
#             for k in range(int(self.numberOfPointRings)):
#
#                 shiftedRing = electrodePositionsOneRing + np.tile(segmentOrientation*self.distanceOfRings*k,
#                                                                   (np.shape(electrodePositionsOneRing)[0],1))
#                 # append shifted ring points to the list of coordinates for this pole
#                 electrodePositionsOneElectrode = np.vstack([electrodePositionsOneElectrode, shiftedRing])
#
#             # if the electrodes are bipolar
#             if self.numberOfPoles == 2:
#                 electrodePositionsPole2 = electrodePositionsOneElectrode + np.tile(segmentOrientation*self.poleDistance, (np.shape(electrodePositionsOneElectrode)[0],1))
#                 electrodePositionsOneElectrode = np.vstack([electrodePositionsOneElectrode, electrodePositionsPole2])
#
#             electrodePositions = np.vstack([electrodePositions, electrodePositionsOneElectrode])
#
#
#         electrodeParameters = {
#                     'sigma' : self.sigma,
#                     'x' : electrodePositions[:,0],  #Coordinates of electrode contacts
#                     'y' : electrodePositions[:,1],
#                     'z' : electrodePositions[:,2],
#                     # 'n' : 20,
#                     # 'r' : 10,
#                     # 'N' : N,
#                     'method': self.method, #or "linesource"
#                 }
#
#         self.electrodeParameters = electrodeParameters
#
# class RecBipolarPoint(RecordingMechanism):
#
#     def __init__(self, radius, numberOfElectrodes=1, poleDistance=5000, positionMax=1, sigma=0.3, method='pointsource'):
#
#         self.radius = radius  # radius of the electrode ring (um)
#         self.poleDistance = poleDistance  # distance between two poles of the same bipolar electrode (um)
#         self.positionMax = positionMax
#
#         super(RecBipolarPoint, self).__init__(numberOfPoints=1, numberOfPoles=2, numberOfElectrodes=numberOfElectrodes, method=method, sigma=sigma)
#
#     def setup_recording_elec(self, bundleGuide, bundleLength):
#
#         bundleGuide = bundleGuide[:, 0:3]
#
#         # first find the bundle guide segment index that corresponds to the intended bundle length (overlap for
#         # myelinated axons gives longer bundle than specified by user
#         bundleLengthIndex = np.shape(bundleGuide)[0] - 1
#         bundleLengthTemp = length_from_coords(bundleGuide)
#         while bundleLengthTemp > bundleLength:
#             bundleLengthIndex -= 1
#             bundleLengthTemp = length_from_coords(bundleGuide[:bundleLengthIndex])
#
#         # lastRecordedSegmentIndex = (np.shape(bundleGuide)[0]-1)*self.positionMax
#         lastRecordedSegmentIndex = bundleLengthIndex * self.positionMax
#
#         if self.numberOfElectrodes > 1:
#             segmentIndices = np.linspace(lastRecordedSegmentIndex / self.numberOfElectrodes,
#                                          lastRecordedSegmentIndex, self.numberOfElectrodes)
#         else:
#             segmentIndices = [lastRecordedSegmentIndex]
#
#         for i in range(self.numberOfElectrodes):
#             self.electrodeDistances.append(np.floor(length_from_coords(bundleGuide[:segmentIndices[i]])))
#
#         electrodePositions = np.array([]).reshape(0, 3)
#
#         for i in range(self.numberOfElectrodes):
#             segmentNumber = int(segmentIndices[i])
#             # segmentNumber = floor(np.shape(bundleGuide)[0]/numberOfElectrodes)*(i+1) - 1
#
#             segmentStartingPos = bundleGuide[segmentNumber - 1, :]
#             segmentEndPos = bundleGuide[segmentNumber, :]
#
#             segmentMiddle = (segmentStartingPos + segmentEndPos) / 2
#             segmentOrientation = segmentEndPos - segmentStartingPos
#             segmentOrientation = segmentOrientation / np.linalg.norm(segmentOrientation)
#
#             # get one random orthogonal vector
#             orthogonalVector = random_perpendicular_vectors(segmentOrientation)[0, :]
#
#             electrodePositionsElectrode1 = (orthogonalVector * self.radius) + \
#                                              (segmentMiddle-segmentOrientation*self.poleDistance/2)
#             electrodePositionsElectrode2 = (orthogonalVector * self.radius) + \
#                                              (segmentMiddle + segmentOrientation * self.poleDistance / 2)
#
#             electrodePositions = np.vstack([electrodePositions, electrodePositionsElectrode1, electrodePositionsElectrode2])
#
#         electrodeParameters = {
#             'sigma': self.sigma,
#             'x': electrodePositions[:, 0],  # Coordinates of electrode contacts
#             'y': electrodePositions[:, 1],
#             'z': electrodePositions[:, 2],
#             # 'n' : 20,
#             # 'r' : 10,
#             # 'N' : N,
#             'method': self.method,  # or "linesource"
#         }
#
#         self.electrodeParameters = electrodeParameters

# TODO: Enable uni-/bipolar recording based on coordinate values for electrode positions
# class RecCoordinates(RecordingMechanism):
#
    # def __init__(self, radius, numberOfElectrodes=1, poleDistance=5000, positionMax=1, sigma=0.3, method='pointsource'):
    #
    #     self.radius = radius  # radius of the electrode ring (um)
    #     self.poleDistance = poleDistance  # distance between two poles of the same bipolar electrode (um)
    #     self.positionMax = positionMax
    #
    #     super(RecCoordinates, self).__init__(numberOfPoints=1, numberOfPoles=2, numberOfElectrodes=numberOfElectrodes,
    #                                           method=method, sigma=sigma)
    #
    # def setup_recording_elec(self, bundleGuide, bundleLength):
    #
    #     bundleGuide = bundleGuide[:, 0:3]
    #
    #     # first find the bundle guide segment index that corresponds to the intended bundle length (overlap for
    #     # myelinated axons gives longer bundle than specified by user
    #     bundleLengthIndex = np.shape(bundleGuide)[0] - 1
    #     bundleLengthTemp = length_from_coords(bundleGuide)
    #     while bundleLengthTemp > bundleLength:
    #         bundleLengthIndex -= 1
    #         bundleLengthTemp = length_from_coords(bundleGuide[:bundleLengthIndex])
    #
    #     # lastRecordedSegmentIndex = (np.shape(bundleGuide)[0]-1)*self.positionMax
    #     lastRecordedSegmentIndex = bundleLengthIndex * self.positionMax
    #
    #     if self.numberOfElectrodes > 1:
    #         segmentIndices = np.linspace(lastRecordedSegmentIndex / self.numberOfElectrodes,
    #                                      lastRecordedSegmentIndex, self.numberOfElectrodes)
    #     else:
    #         segmentIndices = [lastRecordedSegmentIndex]
    #
    #     for i in range(self.numberOfElectrodes):
    #         self.electrodeDistances.append(np.floor(length_from_coords(bundleGuide[:segmentIndices[i]])))
    #
    #     electrodePositions = np.array([]).reshape(0, 3)
    #
    #     for i in range(self.numberOfElectrodes):
    #         segmentNumber = int(segmentIndices[i])
    #         # segmentNumber = floor(np.shape(bundleGuide)[0]/numberOfElectrodes)*(i+1) - 1
    #
    #         segmentStartingPos = bundleGuide[segmentNumber - 1, :]
    #         segmentEndPos = bundleGuide[segmentNumber, :]
    #
    #         segmentMiddle = (segmentStartingPos + segmentEndPos) / 2
    #         segmentOrientation = segmentEndPos - segmentStartingPos
    #         segmentOrientation = segmentOrientation / np.linalg.norm(segmentOrientation)
    #
    #         # get one random orthogonal vector
    #         orthogonalVector = random_perpendicular_vectors(segmentOrientation)[0, :]
    #
    #         electrodePositionsElectrode1 = (orthogonalVector * self.radius) + \
    #                                        (segmentMiddle - segmentOrientation * self.poleDistance / 2)
    #         electrodePositionsElectrode2 = (orthogonalVector * self.radius) + \
    #                                        (segmentMiddle + segmentOrientation * self.poleDistance / 2)
    #
    #         electrodePositions = np.vstack([electrodePositionsElectrode1, electrodePositionsElectrode2])
    #
    #     electrodeParameters = {
    #         'sigma': self.sigma,
    #         'x': electrodePositions[:, 0],  # Coordinates of electrode contacts
    #         'y': electrodePositions[:, 1],
    #         'z': electrodePositions[:, 2],
    #         # 'n' : 20,
    #         # 'r' : 10,
    #         # 'N' : N,
    #         'method': self.method,  # or "linesource"
    #     }
    #
    #     self.electrodeParameters = electrodeParameters


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