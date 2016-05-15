import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d

from axonClass import *
from nameSetters import *

def geometry(bundle):

    if len(bundle.axons[0].xmid) == 0:
        print 'Bundle has not been run yet. No geometry information was generated in NEURON.'
        return

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    axonID = 0
    for axon in bundle.axons:
        if type(axon) == Myelinated:
            style = '--'
        else:
            style = '-'
        ax.plot(axon.xmid, axon.ymid, axon.zmid, style, color= tuple(bundle.axonColors[axonID,:])) # label='axon '+str(axonID),
        ax.text(axon.xmid[-1], axon.ymid[-1], axon.zmid[-1], str(axonID))
        axonID += 1
    ax.plot(bundle.bundleCoords[:,0], bundle.bundleCoords[:,1], bundle.bundleCoords[:,2])#, label='bundle guide')
    #plt.legend()
    plt.title('Axons in space')

    elecCoords = bundle.electrodeCoords
    ax.scatter(elecCoords[:,0], elecCoords[:,1], elecCoords[:,2])

    elecPoles = len(bundle.recordingElecPos)
    for i in range(bundle.numberElecs):
        for j in range(elecPoles):
            # selectionIndices = range(i+j*bundle.numberContactPoints, bundle.numberContactPoints*bundle.numberElecs + j*bundle.numberContactPoints, bundle.numberElecs)
            selectionIndices = range((i*elecPoles+j)*bundle.numberContactPoints, bundle.numberContactPoints*(i*elecPoles+j+1))

            ringCoords = elecCoords[selectionIndices,:]
            ringCoords = np.row_stack((ringCoords, ringCoords[0,:]))
            ax.plot(ringCoords[:,0], ringCoords[:,1], ringCoords[:,2], color=[0.8,0.8,0.8])
            
    
    plt.savefig(os.path.join(bundle.basePath, 'geometry.png'))
    
def CAP1D_singleAxon(bundle, maxNumberOfAxons):

    # get the whole CAP, can be single electrode or multiple
    directory = get_directory_name("CAP1A", bundle.basePath)
    try:
        newestFile = max(glob.iglob(os.path.join(directory,'')+'*.[Dd][Aa][Tt]'), key=os.path.getctime)
    except ValueError:
        print 'No CAP calculation has been performed yet with this set of parameters.'
        return

    CAPraw = np.transpose(np.loadtxt(newestFile))
    time = CAPraw[0,:]
    CAP = CAPraw[1:,:]

    numberOfPlots = min(maxNumberOfAxons, bundle.numberOfAxons)
    axonSelection = np.floor(np.linspace(0,bundle.numberOfAxons-1, numberOfPlots))

    # Subplots
    if len(axonSelection) > 1:
        f, axarr = plt.subplots(numberOfPlots, sharex=True)
    else:
        f = plt.figure()

    for i in range(len(axonSelection)):
        axonIndex = int(axonSelection[i])

        axon = bundle.axons[i]

        axonDiameter = axon.fiberD

        if type(bundle.axons[axonIndex]) == Myelinated:
            axonType = 'myelinated'
        else:
            axonType = 'unmyelinated'

        CAPSingleAxon = CAP[axonIndex,:]

        if len(axonSelection) > 1:
            currentAxis = axarr[i]
        else:
            currentAxis = plt.gca()

        currentAxis.plot(time, CAPSingleAxon, color= tuple(bundle.axonColors[axonIndex,:]))
        currentAxis.set_title('Axon ' + str(axonIndex) + ' (' + axonType + ') with diameter ' + str(axonDiameter) + 'um')
        currentAxis.set_ylabel('CAP [mV]')

    plt.savefig(os.path.join(bundle.basePath, 'CAPSingleAxons.png'))





def CAP1D(bundle, maxNumberOfSubplots = 10):

    # first load the desired data from file
    time, CAP = bundle.get_CAP_from_file()

    numberOfRecordingSites = np.shape(CAP)[0]

    # downsample
    wantedTimeStep = 0.025 #ms
    actualTimeStep = bundle.timeRes
    downsamplingFactor = math.floor(wantedTimeStep/actualTimeStep)
    downsamplingIndexArray = range(0,len(time),int(downsamplingFactor))

    time, CAP = time[downsamplingIndexArray], CAP[:, downsamplingIndexArray]

    if numberOfRecordingSites > 1:

        numberOfPlots = min(maxNumberOfSubplots, numberOfRecordingSites)

        eletrodeSelection = np.floor(np.linspace(0,numberOfRecordingSites-1, numberOfPlots))

        # Subplots
        f, axarr = plt.subplots(numberOfPlots, sharex=True)

        for i in range(numberOfPlots):

            electrodeIndex = eletrodeSelection[i]

            CAPSingleElectrode =  CAP[electrodeIndex,:]
            distanceFromOrigin = bundle.saveParams['L']/numberOfRecordingSites*electrodeIndex

            axarr[i].plot(time, CAPSingleElectrode)
            axarr[i].set_title('distance ' + str(distanceFromOrigin) + ' [um]')
            axarr[i].set_ylabel('CAP [mV]')

            if i == numberOfPlots - 1:
                axarr[i].set_xlabel('time [ms]')
    else:
        fig = plt.figure()
        CAPSingleElectrode =  CAP[numberOfRecordingSites-1,:]
        distanceFromOrigin = bundle.recordingElecPos[0]

        plt.plot(time, CAPSingleElectrode)
        plt.title('distance ' + str(distanceFromOrigin) + ' [um]')
        plt.ylabel('CAP [mV]')
        plt.xlabel('time [ms]')

    plt.savefig(os.path.join(bundle.basePath, 'CAP1D.png'))


def CAP2D(bundle):

    # first load the desired data from file
    time, CAP = bundle.get_CAP_from_file()

    # check if plotting makes sense
    numberOfRecordingSites = np.shape(CAP)[0]

    if numberOfRecordingSites <= 10:
        print 'Plotting of the CAP in two dimensions (time, space) does not make sense with fewer than 10 electrodes. ' \
              'Please select another plotting mechanism or restart the simulation with more electrodes.'
        return


    # print as an image
    fig = plt.figure()
    im = plt.imshow(CAP, cmap=plt.get_cmap('gist_stern'), interpolation='none', aspect='auto')#, norm=LogNorm(vmin=CAPmin, vmax=CAPmax))

    # correct xticks (from samples to ms)
    numberOfXTicks = 10
    tick_locs = np.round(np.linspace(0,np.shape(CAP)[1],numberOfXTicks))
    tick_lbls = np.round(np.linspace(0,bundle.saveParams['tStop'],numberOfXTicks))
    plt.xticks(tick_locs, tick_lbls, fontsize=12)

    # correct yticks (from electrodes to distance)
    numberOfYTicks = 10
    tick_locs = np.round(np.linspace(0,np.shape(CAP)[0],numberOfYTicks))
    tick_lbls = np.round(np.linspace(0,bundle.saveParams['L'],numberOfYTicks))
    plt.yticks(tick_locs, tick_lbls, fontsize=12)

    # add titles, axis labels and colorbar
    fig.suptitle('Compound action potential [uV] over space and time', fontsize=20)
    plt.xlabel('time [ms]')
    plt.ylabel('Distance from axon origin [um]')
    cbar = plt.colorbar(im)

    plt.savefig(os.path.join(bundle.basePath, 'CAP2D.png'))

def voltage(bundle, maxNumberOfSubplots=10):

    timeRec, voltageMatrices = bundle.get_voltage_from_file()

    # now plot
    numberOfAxons = np.shape(voltageMatrices)[0]
    numberOfPlots = min(maxNumberOfSubplots, numberOfAxons)

    axonSelection = np.floor(np.linspace(0,numberOfAxons-1, numberOfPlots))

    if len(axonSelection) > 1:
        f, axarr = plt.subplots(numberOfPlots, sharex=True)
    else:
        f = plt.figure()

    # colors
    jet = plt.get_cmap('jet')

    # for colorbar first check what is the longest axon
    maxAxonLength = 0
    for i in range(len(axonSelection)):
        axonIndex = int(axonSelection[i])
        maxAxonLength = max(maxAxonLength, bundle.axons[axonIndex].L)

    cNorm = colors.Normalize(vmin=0, vmax=int(maxAxonLength)-1)
    scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

    for i in range(len(axonSelection)):
        axonIndex = int(axonSelection[i])

        voltageMatrix = np.transpose(voltageMatrices[axonIndex])

        # find out whether axon is myelinated or not
        isMyelinated = (type(bundle.axons[axonIndex]) == Myelinated)

        axonDiameter = bundle.axons[axonIndex].fiberD
        currentNumberOfSegments = np.shape(voltageMatrix)[1]
        currentAxonLength = bundle.axons[axonIndex].L

        if len(axonSelection) > 1:
            currentAxis = axarr[i]
        else:
            currentAxis = plt.gca()

        if not isMyelinated:
            for j in range(currentNumberOfSegments):
                colorVal = scalarMap.to_rgba(int(j * currentAxonLength / currentNumberOfSegments))
                currentAxis.plot(timeRec, voltageMatrix[:,j], color=colorVal)

            currentAxis.set_title('Voltage of unmyelinated axon with diameter ' + str(axonDiameter) + ' um')
        else:
            Nnodes = bundle.axons[axonIndex].axonnodes

            numberOfRecordingSites = np.shape(voltageMatrix)[1]

            nodePositions = range(0,(Nnodes-1)*11,11)

            nodeDistance = bundle.axons[axonIndex].lengthOneCycle

            nodeCounter = 0
            for j in nodePositions:
                colorVal = scalarMap.to_rgba(int(nodeCounter * nodeDistance))
                currentAxis.plot(np.array(timeRec), np.array(voltageMatrix[:,j]), color=colorVal)
                nodeCounter += 1

            currentAxis.set_title('Voltage of nodes of myelinated axon with diameter ' + str(axonDiameter) + ' um')

        currentAxis.set_ylabel('Voltage [mV]')
        # axarr[i].set_ylim([-100,100])
        currentAxis.set_xlabel('time [ms]')


        # make room for colorbar
        f.subplots_adjust(right=0.8)

        # add colorbar axis
        axColorbar = f.add_axes([0.85, 0.15, 0.05, 0.7])

        cb1 = mpl.colorbar.ColorbarBase(axColorbar, cmap=jet,
                                        norm=cNorm,
                                        orientation='vertical')
        cb1.set_label('length [um]')

    plt.savefig(os.path.join(bundle.basePath, 'voltage.png'))

def voltage_one_myelinated_axon(bundle, myelinatedIndex=0):

    timeRec, voltageMatrices = bundle.get_voltage_from_file()

    # now plot
    numberOfAxons = np.shape(voltageMatrices)[0]
    if myelinatedIndex >= numberOfAxons:
        print 'voltage_one_myelinated_axon: Axon index too high'
        return
    elif myelinatedIndex < 0:
        print 'voltage_one_myelinated_axon: Axon index <0 given.'
        return

    # first search through neurons. This function is for myelinated ones only.
    myelinatedIndexTemp = -1
    axonIndex = 0
    selectedAxon = []
    for axon in bundle.axons:
        # find out whether axon is myelinated or not
        isMyelinated = (type(bundle.axons[axonIndex]) == Myelinated)

        if isMyelinated:
            myelinatedIndexTemp += 1

        if myelinatedIndexTemp == myelinatedIndex:
            selectedAxon = axon
            break

        axonIndex += 1
    if selectedAxon == []:
        print 'Select axon index between 0 and ' + str(myelinatedIndex)
        return


    # extract appropriate voltages from matrix
    voltageMatrix = np.transpose(voltageMatrices[axonIndex])

    # get some fiber properties
    axonDiameter = selectedAxon.fiberD
    # numberOfSegments = np.shape(voltageMatrix)[1]
    axonLength = selectedAxon.L
    Nnodes = bundle.axons[axonIndex].axonnodes
    # numberOfRecordingSites = np.shape(voltageMatrix)[1]
    nodeDistance = bundle.axons[axonIndex].lengthOneCycle

    f, axarr = plt.subplots(11,1)

    # colors
    jet = plt.get_cmap('jet')
    cNorm = colors.Normalize(vmin=0, vmax=int(axonLength)-1)
    scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)


    nodePositions = range(0,(Nnodes-1)*11,11)

    segmentNames = ['node', 'MYSA', 'FLUT', 'STIN', 'STIN', 'STIN', 'STIN', 'STIN', 'STIN', 'FLUT', 'MYSA']

    for i in range(11):

        currentAxis = axarr[i]

        stepFromNode = i

        positions = np.add(nodePositions, stepFromNode) # np.sort(np.concatenate((np.add(nodePositions, 1), np.add(nodePositions, -1))))

        counter = 0
        for j in positions:
            colorVal = scalarMap.to_rgba(int(counter * nodeDistance))
            currentAxis.plot(np.array(timeRec), np.array(voltageMatrix[:,j]), color=colorVal)
            counter += 1

        # currentAxis.set_title('Voltage of segment ' + str(stepFromNode) + ' step(s) behind node of myelinated axon with diameter ' + str(axonDiameter) + ' um')
        if i == 0:
            currentAxis.set_title(segmentNames[i] + ', axon diameter = ' + str(axonDiameter) + 'um')
        else:
            currentAxis.set_title(segmentNames[i])

        currentAxis.set_ylabel('Voltage [mV]')
        currentAxis.set_xlabel('time [ms]')

    # make room for colorbar
    f.subplots_adjust(right=0.8)

    # add colorbar axis
    axColorbar = f.add_axes([0.85, 0.15, 0.05, 0.7])

    cb1 = mpl.colorbar.ColorbarBase(axColorbar, cmap=jet,
                                    norm=cNorm,
                                    orientation='vertical')
    cb1.set_label('length [um]')

    plt.savefig(os.path.join(bundle.basePath, 'voltage.png'))

def diameterHistogram(bundle):
    diametersM = []
    diametersUm = []
    for axon in bundle.axons:
        diameter = axon.fiberD
        if isinstance(axon, Myelinated):
            diametersM.append(diameter)
        else:
            diametersUm.append(diameter)

    # if only unmyelinated axons are present
    if not diametersM and diametersUm:
        plt.hist(diametersUm)
        plt.title('Unmyelinated Diameters')
    # if only myelinated axons are present
    elif diametersM and not diametersUm:
        plt.hist(diametersM)
        plt.title('Myelinated Diameters')
    # if both kinds are present
    elif diametersM and diametersUm:
        f, (ax1, ax2) = plt.subplots(2,1)

        ax1.hist(diametersM)
        ax1.set_title('Myelinated Diameters')

        ax2.hist(diametersUm)
        ax2.set_title('Unmyelinated Diameters')
    # whoops, no axons?
    else:
        print 'No axons in bundle. Cannot plot their diameters.'