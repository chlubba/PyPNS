import os
import glob
import cPickle as pickle

def get_bundle_directory(paramDict, new = False, createDir=False): #dt=0, tStop = 0, pMyel=0, myelinatedDiam = 0, unmyelinatedDiam = 0, L=0, new = False):

    # define here the root of the PyPN file system
    homeDirectory= '/media/carl/4ECC-1C44/PyPN/' # '/home/carl/PNPy/Results/' #  ""#"results"#

    # read out dictionary of parameters (more elegant methon possible?)
    # elecCount = len(paramDict['recordingElecPos'])
    dt=paramDict['timeRes']
    tStop = paramDict['tStop']
    pMyel = paramDict['pMyel']
    L = paramDict['length']
    numberOfAxons = paramDict['numberOfAxons']

    paramsMyel = paramDict['paramsMyel']
    paramsUnmyel = paramDict['paramsUnmyel']

    myelinatedDiam = paramsMyel['fiberD']
    unmyelinatedDiam = paramsUnmyel['fiberD']

    # further process save parameters
    pUnmyel = 1 - pMyel

    # if stimType in ["EXTRA", "INTRA", "NONE"]:
    #     stimulusPathString = "stimType="+stimType+" stimWaveform="+stimWaveform+" stimDutyCycle="+str(stimDutyCycle)+" stimAmplitude="+str(stimAmplitude)+"/"

    if type(myelinatedDiam) in [int, float]:
        myelDiamStr = str(myelinatedDiam)
    else:
        myelDiamStr = 'draw'

    if type(unmyelinatedDiam) in [int, float]:
        unmyelDiamStr = str(unmyelinatedDiam)
    else:
        unmyelDiamStr = 'draw'

    # if elecCount == 1:
    #     poleString = 'monopolarRecording'
    # elif elecCount == 2:
    #     poleString = 'bipolarRecording'
    # else:
    #     print 'Received ' + str(elecCount) + ' as number of poles. Values 1 or 2 allowed only.'
    #     quit()

    #concatenate strings
    pathStringNoStim = "dt="+str(dt)+" tStop="+str(tStop)+" pMyel="+str(pMyel)+" pUnmyel="+str(pUnmyel)+" L="+str(L)+' nAxons='+str(numberOfAxons)#+' '+poleString
    pathString = os.path.join(homeDirectory, pathStringNoStim) # +stimulusPathString

    # find bundle index
    if new:
        versionIndex = 0
        folderName = 'bundle'+str(versionIndex).zfill(5)
        while os.path.isdir(os.path.join(pathString, folderName)):
            versionIndex += 1
            folderName = 'bundle'+str(versionIndex).zfill(5)
        finalBasePath = os.path.join(pathString, folderName)
        os.makedirs(finalBasePath)
    else:
        try:
            latestFolder = max(glob.iglob(os.path.join(pathString,'')+'bundle*'), key=os.path.getmtime)
        except:
            print 'Desired folder empty.'
            return ''

        finalBasePath = latestFolder


    return finalBasePath



def get_directory_name(keyword, basePath):

    """
    retrieve directory name based on the purpose defined by keyword and the bundleParameters

    Args:
        keyword: prefix of the file name
        basePath: bundle basepath

    Returns:

    """


    if keyword=='bundle':
        suffix=''
    else:
        suffix=keyword

    finalCombinedPath = os.path.join(basePath, suffix)

    if not os.path.exists(finalCombinedPath):
            os.makedirs(finalCombinedPath)

    return finalCombinedPath

def get_file_name(recordingType, basePath, newFile=True, directoryType=False):

        if isinstance(directoryType, bool):
            directory = get_directory_name(recordingType, basePath)
        else:
            directory = get_directory_name(directoryType, basePath)

        # filename = 'recording.dat'
        filename = recordingType+'.dat'

        number = 0
        filenameTemp = filename
        if newFile:
            while os.path.isfile(os.path.join(directory,filenameTemp)):
                number += 1
                # print "Be careful this file name already exist! We concatenate a number to the name to avoid erasing your previous file."
                filenameTemp = str(number).zfill(5) + filename


        filename = os.path.join(directory,filenameTemp)


        return filename

def save_bundle(bundle):

    # save the bundle definition file
    bundleSaveLocation = bundle.basePath
    pickle.dump(bundle,open( os.path.join(bundleSaveLocation, 'bundle.cl'), "wb" ))

def open_recent_bundle(Parameters):
    bundleSaveLocation = get_bundle_directory(Parameters, new=False)
    try:
        bundle = pickle.load(open(os.path.join(bundleSaveLocation, 'bundle.cl'), "rb" ))
    except:
        print 'No bundle with these parameters has been generated yet. Set calculationFlag to True.'

    return bundle

def open_bundle_from_location(bundleSaveLocation):

    try:
        bundle = pickle.load(open(os.path.join(bundleSaveLocation, 'bundle.cl'), "rb" ))
    except:
        print 'No bundle with these parameters has been generated yet. Set calculationFlag to True.'

    return bundle