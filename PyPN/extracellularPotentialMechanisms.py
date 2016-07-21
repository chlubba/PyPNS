from abc import abstractmethod, ABCMeta
import numpy as np
from scipy import ndimage
import os
import silencer
import LFPy
import time


class ExtracellularPotentialMechanism(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def calculate_LFP(self, axon, electrodePositions):
        pass

class precomputedFEM(ExtracellularPotentialMechanism):

    def __init__(self, bundleGuide, fieldName='noCuff1'):

        print '\nWhen using a recording mechanism based on a precomputed FEM model, make sure the electrodes are on a ' \
              'long (~1cm) straight part of the bundle guide.'

        fieldDictArray = np.load(os.path.join('/media/carl/4ECC-1C44/ComsolData/usedFields', fieldName, 'fieldDict.npy'))
        self.FEMFieldDict = fieldDictArray[()]

        self.bundleGuide = bundleGuide

    def calculate_LFP(self, axon, electrodePositions):
        """
        Calculate extracellular potential for every electrode point defined in electrodeParameters

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

        def _compute_relative_positions_and_interpolate(axon, electrodePositions):

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
            elecPotentials = np.zeros((electrodePositions.shape[0], np.shape(axon.imem)[1]))

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
                bundleSegStartPointTiled = np.tile(bundleSegStartPoint, (electrodePositions.shape[0], 1))
                elecPosPar = np.inner(electrodePositions - bundleSegStartPointTiled, dir0norm)

                # electrode shifted to origin
                electrodeVector = electrodePositions - (bundleSegStartPointTiled + np.tile(dir0norm,
                                    (electrodePositions.shape[0], 1)) * elecPosPar[:, np.newaxis])

                for axonSegIndLoc in range(len(axonDir2Ds)):
                    axonDir2D = axonDir2Ds[axonSegIndLoc]

                    if isinstance(axonDir2D, int): # if the axon segment lies on the bundle middle exactly
                        elecX = np.ones(electrodePositions.shape[0])*np.linalg.norm(electrodeVector[0,:])
                        elecY = np.zeros(electrodePositions.shape[0])
                    else:
                        axonDir2DTiled = np.tile(axonDir2D, (electrodeVector.shape[0], 1))

                        # electrode coordinates projected onto new base vectors
                        elecX = np.inner(electrodeVector, axonDir2D)

                        elecYVec = electrodeVector - axonDir2DTiled * elecX[:, np.newaxis]
                        elecY = np.linalg.norm(elecYVec, axis=1)

                    axonZ = axonZs[axonSegIndLoc]
                    elecZ = np.add(elecPosPar, -axonZ)

                    axonXDist = np.tile(axonXDists[axonSegIndLoc],(1, len(elecZ)))

                    interpolationPoints = np.vstack([elecX, elecY, elecZ, axonXDist])
                    interpolationPoints = np.divide(interpolationPoints,1000000)  # from um to m TODO: numercal problem?

                    # now interpolate from fieldImage
                    elecPotTempStatic = _interpolateFromImage(self.FEMFieldDict, interpolationPoints, order=1)

                    imemAxonSegInd = axonSegInd-len(axonDir2Ds)+axonSegIndLoc

                    elecPotTemp = np.outer(elecPotTempStatic, axon.imem[imemAxonSegInd,:])*1000 # COMSOL gave V, we need mV


                    # add contributions
                    elecPotentials = np.add(elecPotentials, elecPotTemp)

                bundleSegInd += 1

            return elecPotentials


        # calculate LFP from membrane currents
        return _compute_relative_positions_and_interpolate(axon, electrodePositions)

class homogeneous(ExtracellularPotentialMechanism):

    def __init__(self, method='pointsource', sigma=0.3):

        self.sigma = sigma
        self.method = method


    def calculate_LFP(self, axon, electrodePositions):
        """
        Calculate extracellular potential (LFP) for every electrode point defined in self.electrodeParameters
        Args:
            axon:

        Returns: LFP

        """

        #calculate LFP with LFPy from membrane currents

        # put the locations of electrodes, the method and sigma in the right form for LFPy
        electrodeParameters = {
            'sigma': self.sigma,
            'x': electrodePositions[:, 0],  # Coordinates of electrode contacts
            'y': electrodePositions[:, 1],
            'z': electrodePositions[:, 2],
            'method': self.method,  # 'pointsource' or "linesource"
        }

        # shut down the output, always errors at the end because membrane current too high
        with silencer.nostdout():
            electrodes = LFPy.recextelectrode.RecExtElectrode(axon, **electrodeParameters)

            # calculate LFP by LFPy from membrane current
            electrodes.calc_lfp()


        # get the local field potential
        return electrodes.LFP






