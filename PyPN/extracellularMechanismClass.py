from abc import abstractmethod, ABCMeta
import numpy as np
import os
import silencer
# import LFPy
import time
from extracellularBackend import *

class ExtracellularPotentialMechanism(object):
    __metaclass__ = ABCMeta

    @abstractmethod

    def calculate_extracellular_potential(self, sourcePositions, sourceCurrents, receiverPositions):
        pass

class interpolationExtracellular(ExtracellularPotentialMechanism):

    def __init__(self, bundleGuide, zInterpolator):
        """Stores a precomputed voltage field and calculates the extracelluar potential caused by point sources. Used by :class:`recodingMechanismClass.RecordingMechanism`.

        :param bundleGuide: 3D nerve trajectory
        :param fieldName: string containing the name of the field to import; the location of field files needs to be specified somewhere (TODO)
        """


        self.bundleGuide = bundleGuide
        self.zInterpolator = zInterpolator


    def calculate_extracellular_potential(self, sourcePositions, sourceCurrents, receiverPositions):
        """

        :param sourcePositions: positions of current point sources
        :param sourceCurrents: currents of these point sources
        :param receiverPositions: positions the voltage is calculated for

        """

        def _interpolateStaticPotential(points):
            """
            Here we simplify.

            x = 220 um
            y = 0
            z = 0

            a_x = 0
            a_y = 0
            a_z = only variable

            Args:
                points:
                order:

            Returns:

            """

            a_z = points[-1,:]

            # smoothedFunctionArrayPacked = np.load(os.path.join('/media/carl/4ECC-1C44/ComsolData/interpolated_function_electrode_fixed_a_z_variable', 'smoothed1000.npy'))
            # smoothedFunctionArray1 = smoothedFunctionArrayPacked[()]
            #
            # smoothedFunctionArrayPacked = np.load(
            #     os.path.join('/media/carl/4ECC-1C44/ComsolData/interpolated_function_electrode_fixed_a_z_variable',
            #                  'smoothed.npy'))
            # smoothedFunctionArray0 = smoothedFunctionArrayPacked[()]
            #
            # smoothedFunctionArrayPacked = np.load(
            #     os.path.join('/media/carl/4ECC-1C44/ComsolData/interpolated_function_electrode_fixed_a_z_variable',
            #                  'smoothed1000_zeroed.npy'))
            # smoothedFunctionArray2 = smoothedFunctionArrayPacked[()]
            #
            # smoothedFunctionArrayPacked = np.load(
            #     os.path.join('/media/carl/4ECC-1C44/ComsolData/interpolated_function_electrode_fixed_a_z_variable',
            #                  'smoothed100_zeroed.npy'))
            # smoothedFunctionArray3 = smoothedFunctionArrayPacked[()]

            # import matplotlib.pyplot as plt
            # plt.close('all')
            # plt.figure()
            # plt.plot(smoothedFunctionArray1[0], smoothedFunctionArray1[1])
            # plt.plot(smoothedFunctionArray0[0], smoothedFunctionArray0[1])
            # plt.plot(smoothedFunctionArray1[0], np.maximum(0, np.multiply(np.subtract(1, np.abs(smoothedFunctionArray1[0] / 0.01)), np.max(smoothedFunctionArray0[1]))))
            # plt.show()

            # smoothedInterp = interp1d(smoothedFunctionArray2[0], smoothedFunctionArray3[1], bounds_error=False, fill_value=0)
            #
            # smoothedVoltageStatic = smoothedInterp(a_z)

            # triangularVoltage = np.maximum(0, np.multiply(np.subtract(1, np.abs(a_z / 0.0025)), np.max(smoothedFunctionArray0[1])))


            return self.zInterpolator(a_z) # triangularVoltage # smoothedVoltageStatic #



        def compute_relative_positions_and_interpolate(sourcePositions, sourceCurrents,
                                                                               receiverPositions, currentUnitFEM=-9, currentUnitSource=-9):

            """

            Structure as follows

            for bs in bundleSegments:
                for s in sources on this bundle segment:
                    calculate distance of s from axon guide and distance along axon guide
                for r in all receivers calculate:
                    for s in sources on this bundle segment:
                          calculate receiver position and interpolate potential based on spatial relation between source
                          and receiver and source and bundle

            Args:
                axon: needed for segment positions and membrane currents

            Returns:

            """

            nSourcePoints = np.shape(sourcePositions)[0]

            bundleCoords = self.bundleGuide[:, 0:3]

            # find the bundle guide segments that are closest to the electrode points
            # first calculate bundle segment centers
            bundleSegCenters = bundleCoords[1:, :] - bundleCoords[0:-1, :]

            # then calculate the distances between one electrode point and one bundle segment center
            closestSegInds = []
            receiverDistToOrigins = []
            receiverXDists = []
            for receiverPosition in receiverPositions:
                r2 = (bundleSegCenters[:, 0] - receiverPosition[0]) ** 2 + (bundleSegCenters[:, 1] - receiverPosition[
                    1]) ** 2 + (bundleSegCenters[:, 2] - receiverPosition[2]) ** 2
                r = np.sqrt(r2)

                closestSegInds.append(np.argmin(r))

                bundleSegStartPoint = bundleCoords[closestSegInds[-1]]

                # calculate projection and radius
                dir0 = bundleCoords[closestSegInds[-1] + 1] - bundleSegStartPoint
                dir0norm = dir0 / np.linalg.norm(dir0)

                # distance of receiver along bundle guide
                receiverPosPar = np.inner(receiverPosition - bundleCoords[closestSegInds[-1]], dir0norm)
                receiverDistToOrigin = createGeometry.length_from_coords(
                    bundleCoords[:closestSegInds[-1]]) + receiverPosPar

                receiverDistToOrigins.append(receiverDistToOrigin)

                # normal from bundle guide to axon position
                receiverDir2D = receiverPosition - (bundleSegStartPoint + receiverPosPar * dir0norm)
                receiverXDist = np.linalg.norm(receiverDir2D)

                receiverXDists.append(receiverXDist)

                # if not sourceXDist == 0:
                #     sourceDir2DNorm = sourceDir2D / sourceXDist
                # else:
                #     sourceDir2DNorm = -1

            # first go through all bundle segments, find the associated source positions and calculate the needed
            # quantities for all sources
            receiverPotentials = np.zeros((receiverPositions.shape[0], np.shape(sourceCurrents)[1]))

            # # TODO: delete this again, only for testing
            # lastElectrodeSignal = []
            # f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
            # jet = plt.get_cmap('jet')
            # cNorm = colors.Normalize(vmin=0, vmax=int(nSourcePoints) - 1)
            # scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

            t0 = time.time()

            bundleSegInd = 1
            sourceInd = 0
            sourcesFinishedFlag = False
            lengthAlongBundle = 0
            while bundleSegInd < bundleCoords.shape[
                0] and not sourcesFinishedFlag:  # bundleSegInd points to endpoint of current bundle segment

                bundleSegEndPoint = bundleCoords[bundleSegInd, :]
                bundleSegStartPoint = bundleCoords[bundleSegInd - 1, :]

                # find the normal to the terminating surface of the bundle guide segment
                dir0 = bundleSegEndPoint - bundleSegStartPoint
                dir0norm = dir0 / np.linalg.norm(dir0)
                if bundleSegInd < bundleCoords.shape[0] - 1:
                    dir1 = bundleCoords[bundleSegInd + 1, :] - bundleSegEndPoint
                    n = (dir0 + dir1) / 2
                else:
                    n = dir0

                # process axons associated with the current bundle segment
                # and store axon segment properties in these lists
                sourceDir2Ds = []
                sourceXDists = []
                sourceDistToOrigin = []

                sourcePosition = sourcePositions[sourceInd, :]  # why oh why, no do-while
                # print 'bundleSegEndPoint_ ' + str(bundleSegEndPoint)
                while np.inner(sourcePosition - bundleSegEndPoint, n) < 0:  # while in bundle guide segment
                    # print 'sourcePosition ' + str(sourcePosition)

                    sourcePosition = sourcePositions[sourceInd, :]
                    sourcePosPar = np.inner(sourcePosition - bundleSegStartPoint,
                                            dir0norm)  # compontent parallel to bundle direction

                    # normal from bundle guide to axon position
                    sourceDir2D = sourcePosition - (bundleSegStartPoint + sourcePosPar * dir0norm)
                    sourceXDist = np.linalg.norm(sourceDir2D)
                    if not sourceXDist == 0:
                        sourceDir2DNorm = sourceDir2D / sourceXDist
                    else:
                        sourceDir2DNorm = -1

                    # save all computed values of axon into these lists, they are used when iterating over electrodes
                    sourceXDists.append(sourceXDist)
                    sourceDir2Ds.append(sourceDir2DNorm)
                    sourceDistToOrigin.append(sourcePosPar + lengthAlongBundle)

                    sourceInd += 1

                    if sourceInd >= nSourcePoints:
                        sourcesFinishedFlag = True
                        break

                # now process receivers

                for sourceIndInner, sourceDir2D in enumerate(sourceDir2Ds):
                    # TODO: x,y need to be set correctly, this is only a simplified implementation valid because no side
                    # TODO: displacement of axon segments is modeled
                    receiverX = receiverXDists
                    receiverY = np.zeros(np.shape(receiverXDists))
                    receiverZ = np.zeros(np.shape(receiverXDists))
                    sourceXDist = sourceXDists[sourceIndInner]
                    sourceZDist = np.multiply(receiverDistToOrigins, -1) + sourceDistToOrigin[
                        sourceIndInner]

                    # TODO: always add all three source coordinates.
                    interpolationPoints = np.vstack(
                        [receiverX, receiverY, receiverZ, np.tile(sourceXDist, np.shape(receiverY)), np.zeros(np.shape(receiverY)), sourceZDist])
                    interpolationPoints = np.divide(interpolationPoints,
                                                    1000000)  # from um to m TODO: numercal problem?

                    # now interpolate from fieldImage
                    receiverPotTempStatic = _interpolateStaticPotential(interpolationPoints)

                    imemAxonSegInd = sourceInd - len(sourceDir2Ds) + sourceIndInner

                    # scale potential-voltage-relation with current to obtain temporal signal
                    # COMSOL gave V, we need mV, therefore multiply with 1000
                    # also there can be a mismatch in current unit of the source, eliminate
                    receiverPotTemp = np.outer(receiverPotTempStatic,
                                               sourceCurrents[imemAxonSegInd, :]
                                               * (10) ** (currentUnitSource - currentUnitFEM)) * 1000


                    # add contributions
                    receiverPotentials = np.add(receiverPotentials, receiverPotTemp)

                # measure length along bundle to have the absolute position along it for all segments.
                lengthAlongBundle += np.linalg.norm(bundleSegEndPoint - bundleSegStartPoint)

                # go to next segment
                bundleSegInd += 1

            return receiverPotentials

        # calculate LFP from membrane currents
        return compute_relative_positions_and_interpolate(sourcePositions, sourceCurrents, receiverPositions)

class precomputedFEM_6D(ExtracellularPotentialMechanism):

    def __init__(self, bundleGuide, fieldName='oil_different_positions_2cm', receiverDisplacement=0, interpolationOrder=1):
        """Stores a precomputed voltage field and calculates the extracelluar potential caused by point sources. Used by :class:`recodingMechanismClass.RecordingMechanism`.

        :param bundleGuide: 3D nerve trajectory
        :param fieldName: string containing the name of the field to import; the location of field files needs to be specified somewhere (TODO)
        """

        # TODO: finally solve the location issue, saving and field loading. Input dict to bundle.

        # todo: input params that allow more degrees of freedom for the FEM model

        print '\nWhen using a recording mechanism based on a precomputed FEM model, make sure the electrodes are on a ' \
              'long (~1cm) straight part of the bundle guide.'

        fieldDictArray = np.load(
            os.path.join('/media/carl/4ECC-1C44/ComsolData/usedFields', fieldName, 'fieldDictOrdered.npy'))
        self.FEMFieldDict = fieldDictArray[()]

        self.bundleGuide = bundleGuide
        self.receiverDisplacement = receiverDisplacement
        self.order = interpolationOrder

        # parse for symmetries in the FEM simulation. Assume the source positions were chosen to minimize number of
        # calculations
        self.symmetries = np.zeros(3)
        for dimensionIndex, axonKey in enumerate(['axonX', 'axonY', 'axonZ']):

            axonPositionValues = self.FEMFieldDict[axonKey]

            # if only one direction was explored, this shows at least one symmetry axis perpendiacular to the axis
            # extend is present
            if np.min(axonPositionValues) == 0:

                self.symmetries[dimensionIndex] = 1

                # if only one position is exported, two symmetry axes seem to exist
                if np.max(axonPositionValues) == 0:

                    self.symmetries[dimensionIndex] = 2


    def calculate_extracellular_potential(self, sourcePositions, sourceCurrents, receiverPositions):
        """

        :param sourcePositions: positions of current point sources
        :param sourceCurrents: currents of these point sources
        :param receiverPositions: positions the voltage is calculated for

        """

        def _getImageCoords1D(physicalCoordinate1D, pointsOfInterest):

            if len(physicalCoordinate1D) > 1:

                # interpolate to obtain coordinate position
                physicalCoordinate1D.sort()
                # if out of bounds, give a random value outside the interpolation range of image to obtain 0 as output
                coordInterp = interp1d(physicalCoordinate1D, range(len(physicalCoordinate1D)),
                                       bounds_error=False, fill_value=physicalCoordinate1D[0]-1)

                pointsOfInterest = np.array(pointsOfInterest)

                coords = coordInterp(pointsOfInterest)

            else:
                # a single value only signifies that this coordinate does not interest. Only for source positions.
                coords = np.zeros(len(pointsOfInterest))

            return coords

        # # todo: uncomment the version underneath. This is only for testing purposes.
        # def _interpolateStaticPotentialFromImage(points, order=1):
        #     """
        #     Here we simplify.
        #
        #     x = 220 um
        #     y = 0
        #     z = 0
        #
        #     a_x = 0
        #     a_y = 0
        #     a_z = only variable
        #
        #     Args:
        #         points:
        #         order:
        #
        #     Returns:
        #
        #     """
        #
        #     a_z = points[-1,:]
        #
        #     smoothedFunctionArrayPacked = np.load(os.path.join('/media/carl/4ECC-1C44/ComsolData/interpolated_function_electrode_fixed_a_z_variable', 'smoothed1000.npy'))
        #     smoothedFunctionArray1 = smoothedFunctionArrayPacked[()]
        #
        #     smoothedFunctionArrayPacked = np.load(
        #         os.path.join('/media/carl/4ECC-1C44/ComsolData/interpolated_function_electrode_fixed_a_z_variable',
        #                      'smoothed.npy'))
        #     smoothedFunctionArray0 = smoothedFunctionArrayPacked[()]
        #
        #     smoothedFunctionArrayPacked = np.load(
        #         os.path.join('/media/carl/4ECC-1C44/ComsolData/interpolated_function_electrode_fixed_a_z_variable',
        #                      'smoothed1000_zeroed.npy'))
        #     smoothedFunctionArray2 = smoothedFunctionArrayPacked[()]
        #
        #     smoothedFunctionArrayPacked = np.load(
        #         os.path.join('/media/carl/4ECC-1C44/ComsolData/interpolated_function_electrode_fixed_a_z_variable',
        #                      'smoothed100_zeroed.npy'))
        #     smoothedFunctionArray3 = smoothedFunctionArrayPacked[()]
        #
        #     # import matplotlib.pyplot as plt
        #     # plt.close('all')
        #     # plt.figure()
        #     # plt.plot(smoothedFunctionArray1[0], smoothedFunctionArray1[1])
        #     # plt.plot(smoothedFunctionArray0[0], smoothedFunctionArray0[1])
        #     # plt.plot(smoothedFunctionArray1[0], np.maximum(0, np.multiply(np.subtract(1, np.abs(smoothedFunctionArray1[0] / 0.01)), np.max(smoothedFunctionArray0[1]))))
        #     # plt.show()
        #
        #     # smoothedInterp = interp1d(smoothedFunctionArray2[0], smoothedFunctionArray3[1], bounds_error=False, fill_value=0)
        #     #
        #     # smoothedVoltageStatic = smoothedInterp(a_z)
        #
        #     triangularVoltage = np.maximum(0, np.multiply(np.subtract(1, np.abs(a_z / 0.0025)), np.max(smoothedFunctionArray0[1])))
        #
        #
        #     return triangularVoltage # smoothedVoltageStatic #

        def _interpolateStaticPotentialFromImage(points, order=1):

            fieldKeys = (('axonX', 'x'), ('axonY', 'y'), ('axonZ', 'z'))
            receiverImageCoordsMat = []
            sourceImageCoordsMat = []

            # process coordinate pairs (axon position and electrode position). If for a coordinate the axon position is only
            # given in one direction (only positive coordinate values), for negative coordinate values the field of the positive
            # position is taken. But then also the receiver position needs to be mirrored on the axis of symmetry.
            for keyIndex, coordKeyPair in enumerate(fieldKeys):

                # get the two coordinate vectors (receiver and source)
                receiverCoords = points[keyIndex,:]
                sourceCoords = points[keyIndex + 3, :]

                # check which symmetries exist and manipulate coordinates accordingly. If no symmetry
                # (self.symmetries[keyIndex] == 0), no action is neccessary.

                # first possibility, one symmetry axis
                if self.symmetries[keyIndex] == 1:

                    # mirror source positions at symmetry axis
                    sourceSign = np.sign(sourceCoords)
                    sourceSign[sourceSign == 0] = 1
                    receiverCoords = np.multiply(receiverCoords, sourceSign)

                    # take absolute value (negative positions have the same field but mirrored, therefore a mirroring
                    # of the receiver positions)
                    sourceCoords = np.abs(sourceCoords)

                elif self.symmetries[keyIndex] == 2:

                    # position changes of the source along this axis can be captured otherwise (in another axis)

                    # as it is irrelevant it will be set to zero
                    sourceCoords = np.zeros(np.shape(sourceCoords))

                    # the field is bound to be symmetry along this axis, so the receiver position will only be taken in
                    # one direction
                    receiverCoords = np.abs(receiverCoords)

                receiverImageCoords = _getImageCoords1D(self.FEMFieldDict[coordKeyPair[1]], receiverCoords)
                sourceImageCoords = _getImageCoords1D(self.FEMFieldDict[coordKeyPair[0]], sourceCoords)

                # append them to their respective matrix to enable the right order
                receiverImageCoordsMat.append(receiverImageCoords)
                sourceImageCoordsMat.append(sourceImageCoords)

            # concatenate matrices
            combinedImageCoords = np.vstack((np.array(receiverImageCoordsMat), np.array(sourceImageCoordsMat)))

            return ndimage.map_coordinates(self.FEMFieldDict['fieldImage'], combinedImageCoords, order=order)

        # todo: only for testing

        import matplotlib.pyplot as plt
        plt.figure()

        # for zTest in [0]: # np.linspace(0, 0.01, 10):
        #
        #     x = [0.00022]
        #     y = [0]
        #     z = [zTest]
        #     axonX = [0]
        #     axonY = [0]
        #     axonZ = np.linspace(-0.04, 0.04, 20000)
        #
        #     testCoords = np.vstack((np.tile(x, np.shape(axonZ)), np.tile(y, np.shape(axonZ)), np.tile(z, np.shape(axonZ)),
        #                             np.tile(axonX, np.shape(axonZ)), np.tile(axonY, np.shape(axonZ)), axonZ))
        #     voltageValues = _interpolateStaticPotentialFromImage(testCoords, order=2)
        #
        #     plt.plot(axonZ, voltageValues, label='second order interpolation')

        def smooth(y, box_pts):
            box = np.ones(box_pts) / box_pts
            y_smooth = np.convolve(y, box, mode='same')
            return y_smooth

        for zTest in [0]: # np.linspace(0, 0.01, 10):

            # x = [0.00022]
            # y = [0]
            # z = [zTest]
            # axonX = [0]
            # axonY = [0]
            # axonZ = np.linspace(-0.04, 0.04, 20000)
            #
            # testCoords = np.vstack((np.tile(x, np.shape(axonZ)), np.tile(y, np.shape(axonZ)), np.tile(z, np.shape(axonZ)),
            #                         np.tile(axonX, np.shape(axonZ)), np.tile(axonY, np.shape(axonZ)), axonZ))
            # voltageValues = _interpolateStaticPotentialFromImage(testCoords, order=1)
            #
            # plt.plot(axonZ, smooth(voltageValues, 100), label='smoothed 100 points')
            # plt.plot(axonZ, smooth(voltageValues, 500), label='smoothed 500 points')
            # plt.plot(axonZ, smooth(voltageValues, 1000), label='smoothed 1000 points')
            #
            # plt.plot(axonZ, smooth(smooth(voltageValues, 10), 10), label='smoothed twice 10 points')
            # plt.plot(axonZ, smooth(smooth(voltageValues, 50), 50), label='smoothed twice 50 points')
            #
            # plt.plot(axonZ, smooth(smooth(smooth(voltageValues, 50), 50), 50), label='smoothed three times 50 points')
            #
            # smoothedCurve = []
            # smoothedCurve.append(axonZ)
            # smoothedCurve.append(smooth(voltageValues, 1000))

            x = [0.00022]
            y = [0]
            z = [zTest]
            axonX = [0]
            axonY = [0]
            axonZ = [-0.0105, -0.01, -0.0075, -0.005, -0.0025, 0, 0.0025, 0.005, 0.0075, 0.01, 0.0105]

            testCoords = np.vstack(
                (np.tile(x, np.shape(axonZ)), np.tile(y, np.shape(axonZ)), np.tile(z, np.shape(axonZ)),
                 np.tile(axonX, np.shape(axonZ)), np.tile(axonY, np.shape(axonZ)), axonZ))
            voltageValues = _interpolateStaticPotentialFromImage(testCoords, order=1)

            axonZ = np.concatenate(([-0.011], axonZ, [0.011]))
            voltageValues = np.concatenate(([0], voltageValues, [0]))

            plt.plot(axonZ, voltageValues, label='linear')

            voltageInterpolator = interp1d(axonZ, voltageValues, bounds_error=False, fill_value=0)
            axonZInterp = np.linspace(-0.02, 0.02, 10000)
            voltageInterp = voltageInterpolator(axonZInterp)


            plt.plot(axonZInterp, smooth(voltageInterp, 100), label='smoothed 100 points')
            plt.plot(axonZInterp, smooth(voltageInterp, 500), label='smoothed 500 points')
            plt.plot(axonZInterp, smooth(voltageInterp, 1000), label='smoothed 1000 points')

            plt.plot(axonZInterp, smooth(smooth(voltageInterp, 10), 10), label='smoothed twice 10 points')
            plt.plot(axonZInterp, smooth(smooth(voltageInterp, 50), 50), label='smoothed twice 50 points')

            plt.plot(axonZInterp, smooth(smooth(smooth(voltageInterp, 50), 50), 50), label='smoothed three times 50 points')

            smoothedCurve = []
            smoothedCurve.append(axonZInterp)
            smoothedCurve.append(voltageInterp)

            np.save(os.path.join('/media/carl/4ECC-1C44/ComsolData/interpolated_function_electrode_fixed_a_z_variable', 'smoothed1_zeroed.npy'), smoothedCurve)

        plt.legend()
        plt.show()

        def compute_relative_positions_and_interpolate_symmetric_inhomogeneity(sourcePositions, sourceCurrents,
                                                                               receiverPositions, currentUnitFEM=-9, currentUnitSource=-9):

            """

            Structure as follows

            for bs in bundleSegments:
                for s in sources on this bundle segment:
                    calculate distance of s from axon guide and distance along axon guide
                for r in all receivers calculate:
                    for s in sources on this bundle segment:
                          calculate receiver position and interpolate potential based on spatial relation between source
                          and receiver and source and bundle

            Args:
                axon: needed for segment positions and membrane currents

            Returns:

            """

            nSourcePoints = np.shape(sourcePositions)[0]

            bundleCoords = self.bundleGuide[:, 0:3]

            # find the bundle guide segments that are closest to the electrode points
            # first calculate bundle segment centers
            bundleSegCenters = bundleCoords[1:, :] - bundleCoords[0:-1, :]

            # then calculate the distances between one electrode point and one bundle segment center
            closestSegInds = []
            receiverDistToOrigins = []
            receiverXDists = []
            for receiverPosition in receiverPositions:
                r2 = (bundleSegCenters[:, 0] - receiverPosition[0]) ** 2 + (bundleSegCenters[:, 1] - receiverPosition[
                    1]) ** 2 + (bundleSegCenters[:, 2] - receiverPosition[2]) ** 2
                r = np.sqrt(r2)

                closestSegInds.append(np.argmin(r))

                bundleSegStartPoint = bundleCoords[closestSegInds[-1]]

                # calculate projection and radius
                dir0 = bundleCoords[closestSegInds[-1] + 1] - bundleSegStartPoint
                dir0norm = dir0 / np.linalg.norm(dir0)

                # distance of receiver along bundle guide
                receiverPosPar = np.inner(receiverPosition - bundleCoords[closestSegInds[-1]], dir0norm)
                receiverDistToOrigin = createGeometry.length_from_coords(
                    bundleCoords[:closestSegInds[-1]]) + receiverPosPar

                receiverDistToOrigins.append(receiverDistToOrigin)

                # normal from bundle guide to axon position
                receiverDir2D = receiverPosition - (bundleSegStartPoint + receiverPosPar * dir0norm)
                receiverXDist = np.linalg.norm(receiverDir2D)

                receiverXDists.append(receiverXDist)

                # if not sourceXDist == 0:
                #     sourceDir2DNorm = sourceDir2D / sourceXDist
                # else:
                #     sourceDir2DNorm = -1

            # first go through all bundle segments, find the associated source positions and calculate the needed
            # quantities for all sources
            receiverPotentials = np.zeros((receiverPositions.shape[0], np.shape(sourceCurrents)[1]))

            # # TODO: delete this again, only for testing
            # lastElectrodeSignal = []
            # f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
            # jet = plt.get_cmap('jet')
            # cNorm = colors.Normalize(vmin=0, vmax=int(nSourcePoints) - 1)
            # scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

            t0 = time.time()

            bundleSegInd = 1
            sourceInd = 0
            sourcesFinishedFlag = False
            lengthAlongBundle = 0
            while bundleSegInd < bundleCoords.shape[
                0] and not sourcesFinishedFlag:  # bundleSegInd points to endpoint of current bundle segment

                bundleSegEndPoint = bundleCoords[bundleSegInd, :]
                bundleSegStartPoint = bundleCoords[bundleSegInd - 1, :]

                # find the normal to the terminating surface of the bundle guide segment
                dir0 = bundleSegEndPoint - bundleSegStartPoint
                dir0norm = dir0 / np.linalg.norm(dir0)
                if bundleSegInd < bundleCoords.shape[0] - 1:
                    dir1 = bundleCoords[bundleSegInd + 1, :] - bundleSegEndPoint
                    n = (dir0 + dir1) / 2
                else:
                    n = dir0

                # process axons associated with the current bundle segment
                # and store axon segment properties in these lists
                sourceDir2Ds = []
                sourceXDists = []
                sourceDistToOrigin = []

                sourcePosition = sourcePositions[sourceInd, :]  # why oh why, no do-while
                # print 'bundleSegEndPoint_ ' + str(bundleSegEndPoint)
                while np.inner(sourcePosition - bundleSegEndPoint, n) < 0:  # while in bundle guide segment
                    # print 'sourcePosition ' + str(sourcePosition)

                    sourcePosition = sourcePositions[sourceInd, :]
                    sourcePosPar = np.inner(sourcePosition - bundleSegStartPoint,
                                            dir0norm)  # compontent parallel to bundle direction

                    # normal from bundle guide to axon position
                    sourceDir2D = sourcePosition - (bundleSegStartPoint + sourcePosPar * dir0norm)
                    sourceXDist = np.linalg.norm(sourceDir2D)
                    if not sourceXDist == 0:
                        sourceDir2DNorm = sourceDir2D / sourceXDist
                    else:
                        sourceDir2DNorm = -1

                    # save all computed values of axon into these lists, they are used when iterating over electrodes
                    sourceXDists.append(sourceXDist)
                    sourceDir2Ds.append(sourceDir2DNorm)
                    sourceDistToOrigin.append(sourcePosPar + lengthAlongBundle)

                    sourceInd += 1

                    if sourceInd >= nSourcePoints:
                        sourcesFinishedFlag = True
                        break

                # now process receivers

                for sourceIndInner, sourceDir2D in enumerate(sourceDir2Ds):
                    # TODO: x,y need to be set correctly, this is only a simplified implementation valid because no side
                    # TODO: displacement of axon segments is modeled
                    receiverX = receiverXDists
                    receiverY = np.zeros(np.shape(receiverXDists))
                    receiverZ = np.zeros(np.shape(receiverXDists))
                    sourceXDist = sourceXDists[sourceIndInner]
                    sourceZDist = np.multiply(receiverDistToOrigins, -1) + sourceDistToOrigin[
                        sourceIndInner] + self.receiverDisplacement

                    # TODO: always add all three source coordinates.
                    interpolationPoints = np.vstack(
                        [receiverX, receiverY, receiverZ, np.tile(sourceXDist, np.shape(receiverY)), np.zeros(np.shape(receiverY)), sourceZDist])
                    interpolationPoints = np.divide(interpolationPoints,
                                                    1000000)  # from um to m TODO: numercal problem?

                    # now interpolate from fieldImage
                    receiverPotTempStatic = _interpolateStaticPotentialFromImage(interpolationPoints, order=self.order)

                    imemAxonSegInd = sourceInd - len(sourceDir2Ds) + sourceIndInner

                    # scale potential-voltage-relation with current to obtain temporal signal
                    # COMSOL gave V, we need mV, therefore multiply with 1000
                    # also there can be a mismatch in current unit of the source, eliminate
                    receiverPotTemp = np.outer(receiverPotTempStatic,
                                               sourceCurrents[imemAxonSegInd, :]
                                               * (10) ** (currentUnitSource - currentUnitFEM)) * 1000

                    # # compontent parallel to bundle direction
                    # bundleSegStartPointTiled = np.tile(bundleSegStartPoint, (receiverPositions.shape[0], 1))
                    # receiverPosPar = np.inner(receiverPositions - bundleSegStartPointTiled, dir0norm)
                    #
                    # # electrode shifted to origin
                    # receiverVector = receiverPositions - \
                    #                  (bundleSegStartPointTiled + np.tile(dir0norm,(receiverPositions.shape[0], 1)) * receiverPosPar[:,np.newaxis])

                    # for sourceLoc in range(len(sourceDir2Ds)):
                    #     sourceDir2D = sourceDir2Ds[sourceLoc]
                    #
                    #     if isinstance(sourceDir2D, int):  # if the axon segment lies on the bundle middle exactly
                    #         receiverX = np.ones(receiverPositions.shape[0]) * np.linalg.norm(receiverVector[0, :])
                    #         receiverY = np.zeros(receiverPositions.shape[0])
                    #     else:
                    #         sourceDir2DTiled = np.tile(sourceDir2D, (receiverVector.shape[0], 1))
                    #
                    #         # electrode coordinates projected onto new base vectors
                    #         receiverX = np.inner(receiverVector, sourceDir2D)
                    #
                    #         receiverYVec = receiverVector - sourceDir2DTiled * receiverX[:, np.newaxis]
                    #         receiverY = np.linalg.norm(receiverYVec, axis=1)
                    #
                    #     sourceDistToOrigin = sourceZs[sourceLoc]
                    #     receiverZ = np.add(receiverPosPar, -sourceZ)
                    #
                    #     sourceXDist = np.tile(sourceXDists[sourceLoc], (1, len(receiverZ)))
                    #
                    #     interpolationPoints = np.vstack([receiverX, receiverY, receiverZ, sourceXDist])
                    #     interpolationPoints = np.divide(interpolationPoints,
                    #                                     1000000)  # from um to m TODO: numercal problem?
                    #
                    #     # now interpolate from fieldImage
                    #     receiverPotTempStatic = _interpolateFromImage(fieldDict, interpolationPoints, order=1)
                    #
                    #     imemAxonSegInd = sourceInd - len(sourceDir2Ds) + sourceLoc
                    #
                    #     # scale potential-voltage-relation with current to obtain temporal signal
                    #     # COMSOL gave V, we need mV, therefore multiply with 1000
                    #     # also there can be a mismatch in current unit of the source, eliminate
                    #     receiverPotTemp = np.outer(receiverPotTempStatic,
                    #                            sourceCurrents[imemAxonSegInd, :]
                    #                            * (10)**(currentUnitSource-currentUnitFEM)) * 1000
                    #
                    #     # import matplotlib.pyplot as plt
                    #     # plt.plot(receiverPotTemp.T)
                    #     # plt.show()
                    #
                    # add contributions
                    receiverPotentials = np.add(receiverPotentials, receiverPotTemp)

                # measure length along bundle to have the absolute position along it for all segments.
                lengthAlongBundle += np.linalg.norm(bundleSegEndPoint - bundleSegStartPoint)

                # go to next segment
                bundleSegInd += 1

            return receiverPotentials

        # calculate LFP from membrane currents
        return compute_relative_positions_and_interpolate_symmetric_inhomogeneity(sourcePositions, sourceCurrents, receiverPositions)



class precomputedFEM_symmetrical_inhomogeneity(ExtracellularPotentialMechanism):

    def __init__(self, bundleGuide, fieldName='oil_different_positions_2cm', receiverDisplacement=0):
        """Stores a precomputed voltage field and calculates the extracelluar potential caused by point sources. Used by :class:`recodingMechanismClass.RecordingMechanism`.

        :param bundleGuide: 3D nerve trajectory
        :param fieldName: string containing the name of the field to import; the location of field files needs to be specified somewhere (TODO)
        """

        # TODO: finally solve the location issue, saving and field loading. Input dict to bundle.

        # todo: input params that allow more degrees of freedom for the FEM model

        print '\nWhen using a recording mechanism based on a precomputed FEM model, make sure the electrodes are on a ' \
              'long (~1cm) straight part of the bundle guide.'

        fieldDictArray = np.load(
            os.path.join('/media/carl/4ECC-1C44/ComsolData/usedFields', fieldName, 'fieldDict.npy'))
        self.FEMFieldDict = fieldDictArray[()]

        self.bundleGuide = bundleGuide
        self.receiverDisplacement = receiverDisplacement

    # def calculate_LFP(self, axon, electrodePositions): # sourcePositions, sourceCurrents,
    def calculate_extracellular_potential(self, sourcePositions, sourceCurrents, receiverPositions):
        """

        :param sourcePositions: positions of current point sources
        :param sourceCurrents: currents of these point sources
        :param receiverPositions: positions the voltage is calculated for

        """

        # calculate LFP from membrane currents
        return compute_relative_positions_and_interpolate_symmetric_inhomogeneity(sourcePositions, sourceCurrents, receiverPositions, self.FEMFieldDict, self.bundleGuide, self.receiverDisplacement)

class precomputedFEM(ExtracellularPotentialMechanism):

    def __init__(self, bundleGuide, fieldName='noCuff1'):
        """Stores a precomputed voltage field and calculates the extracelluar potential caused by point sources. Used by :class:`recodingMechanismClass.RecordingMechanism`.

        :param bundleGuide: 3D nerve trajectory
        :param fieldName: string containing the name of the field to import; the location of field files needs to be specified somewhere (TODO)
        """

        # TODO: finally solve the location issue, saving and field loading. Input dict to bundle.

        # todo: input params that allow more degrees of freedom for the FEM model

        print '\nWhen using a recording mechanism based on a precomputed FEM model, make sure the electrodes are on a ' \
              'long (~1cm) straight part of the bundle guide.'

        fieldDictArray = np.load(
            os.path.join('/media/carl/4ECC-1C44/ComsolData/usedFields', fieldName, 'fieldDict.npy'))
        self.FEMFieldDict = fieldDictArray[()]

        self.bundleGuide = bundleGuide

    # def calculate_LFP(self, axon, electrodePositions): # sourcePositions, sourceCurrents,
    def calculate_extracellular_potential(self, sourcePositions, sourceCurrents, receiverPositions):
        """

        :param sourcePositions: positions of current point sources
        :param sourceCurrents: currents of these point sources
        :param receiverPositions: positions the voltage is calculated for

        """

        # calculate LFP from membrane currents
        return compute_relative_positions_and_interpolate(sourcePositions, sourceCurrents, receiverPositions, self.FEMFieldDict, self.bundleGuide)



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






