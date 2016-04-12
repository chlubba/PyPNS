import os
import glob

def get_bundle_directory(paramDict, new = False): #dt=0, tStop = 0, p_A=0, myelinatedDiam = 0, unmyelinatedDiam = 0, L=0, new = False):

    elecCount = len(paramDict['recordingElecPos'])
    dt=paramDict['timeRes']
    tStop = paramDict['tStop']
    p_A = paramDict['p_A']
    L = paramDict['lengthOfBundle']
    numberOfAxons = paramDict['numberOfAxons']

    myelinated_A = paramDict['myelinated_A']
    unmyelinated = paramDict['unmyelinated']

    unmyelinatedDiam = unmyelinated['fiberD']
    myelinatedDiam = myelinated_A['fiberD']


    homeDirectory="/media/carl/4ECC-1C44/PyPN/"#""#"results"#

    # prepare single parameter values for string insertion
    p_C = 1 - p_A

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

    if elecCount == 1:
        poleString = 'monopolarRecording'
    elif elecCount == 2:
        poleString = 'bipolarRecording'
    else:
        print 'Received ' + str(elecCount) + ' as number of poles. Values 1 or 2 allowed only.'
        quit()

    #concatenate strings
    pathStringNoStim = "dt="+str(dt)+" tStop="+str(tStop)+" p_A="+str(p_A)+" p_C="+str(p_C)+" L="+str(L)+' nAxons='+str(numberOfAxons)+' '+poleString
    pathString = os.path.join(homeDirectory, pathStringNoStim) # +stimulusPathString

    # find bundle index
    if new:
        versionIndex = 0
        folderName = 'bundle'+str(versionIndex).zfill(5)
        while os.path.isdir(os.path.join(pathString, folderName)):
            versionIndex += 1
            folderName = 'bundle'+str(versionIndex).zfill(5)
        finalBasePath = os.path.join(pathString, folderName)
    else:
        try:
            latestFolder = max(glob.iglob(os.path.join(pathString,'')+'bundle*'), key=os.path.getmtime)
        except:
            print 'Desired folder empty.'
            return ''

        finalBasePath = latestFolder

    return finalBasePath



def get_directory_name(keyword, basePath):
    # retrieve directory name based on the purpose defined by keyword and the bundleParameters
    # 4 possible keywords:
    # elec -> temporary electrode folder
    # draw -> distrbution folder
    # CAP -> compound action potential folder
    # V -> voltage folder

    suffix = {
        'elec': "electrodes",
        'draw': "draws",
        'CAP': "CAP",
        'CAP1A': "CAPSingleAxons",
        'V': "Voltage",
        'bundle' : ""
    }.get(keyword,-1)

    finalCombinedPath = os.path.join(basePath, suffix)

    if not os.path.exists(finalCombinedPath):
            os.makedirs(finalCombinedPath)

    return finalCombinedPath

def get_file_name(recordingType, basePath, newFile=True):

        directory = get_directory_name(recordingType, basePath)

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
