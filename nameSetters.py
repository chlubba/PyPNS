import os
import glob

def getBundleDirectory(dt=0, tStop = 0, p_A=0, myelinatedDiam = 0, unmyelinatedDiam = 0, L=0, elecCount=2, stimType = "EXTRA", stimWaveform = "", stimDutyCycle = 0, stimAmplitude = 0, new = False):

    homeDirectory="/media/carl/4ECC-1C44/PyPN/"#""#

    # prepare single parameter values for string insertion
    p_C = 1 - p_A

    if stimType in ["EXTRA", "INTRA", "NONE"]:
        stimulusPathString = "stimType="+stimType+" stimWaveform="+stimWaveform+" stimDutyCycle="+str(stimDutyCycle)+" stimAmplitude="+str(stimAmplitude)+"/"

    if type(myelinatedDiam) in [int, float]:
        myelDiamStr = str(myelinatedDiam)
    else:
        myelDiamStr = 'draw'

    if type(unmyelinatedDiam) in [int, float]:
        unmyelDiamStr = str(unmyelinatedDiam)
    else:
        unmyelDiamStr = 'draw'

    #concatenate strings
    pathStringNoStim = "dt="+str(dt)+" tStop="+str(tStop)+" p_A="+str(p_A)+" p_C="+str(p_C)+" L="+str(L)+"/"#+" myelinatedDiam="+myelDiamStr+" unmyelinatedDiam="+unmyelDiamStr
    pathString = homeDirectory+stimulusPathString+pathStringNoStim

    # find bundle index
    if new:
        versionIndex = 0
        folderName = 'bundle'+str(versionIndex).zfill(5)+'/'
        while os.path.isdir(pathString+folderName):
            versionIndex += 1
            folderName = 'bundle'+str(versionIndex).zfill(5)+'/'
        finalBasePath = pathString+folderName
    else:
        try:
            latestFolder = max(glob.iglob(pathString+'bundle*'), key=os.path.getmtime)
        except:
            print 'Desired folder empty.'
            return ''

        finalBasePath = latestFolder + '/'

    return finalBasePath



def getDirectoryName(keyword, basePath):
    # retrieve directory name based on the purpose defined by keyword and the bundleParameters
    # 4 possible keywords:
    # elec -> temporary electrode folder
    # draw -> distrbution folder
    # CAP -> compound action potential folder
    # V -> voltage folder

    suffix = {
        'elec': "electrodes/",
        'draw': "draws/",
        'CAP': "CAP/",
        'CAP1A': "CAPSingleAxons/",
        'V': "Voltage/",
        'bundle' : ""
    }.get(keyword,-1)

    finalCombinedPath = basePath + suffix

    if not os.path.exists(finalCombinedPath):
            os.makedirs(finalCombinedPath)

    return finalCombinedPath

def getFileName(recordingType, basePath):

        directory = getDirectoryName(recordingType, basePath)

        # filename = 'recording.dat'
        filename = recordingType+'.dat'

        number = 0
        filenameTemp = filename
        while os.path.isfile(directory+filenameTemp):
            number += 1
            # print "Be careful this file name already exist! We concatenate a number to the name to avoid erasing your previous file."
            filenameTemp = str(number).zfill(5) + filename

        filename = directory+filenameTemp

        return filename