import os

def getDirectoryName(keyword, dt=0, tStop = 0, p_A=0, myelinatedDiam = 0, unmyelinatedDiam = 0, L=0, elecCount=2, stimType = "EXTRA", stimWaveform = "", stimDutyCycle = 0, stimAmplitude = 0):
    # retrieve directory name based on the purpose defined by keyword and the bundleParameters
    # 4 possible keywords:
    # elec -> temporary electrode folder
    # draw -> distrbution folder
    # CAP -> compound action potential folder
    # V -> voltage folder

    homeDirectory="/media/carl/4ECC-1C44/PyPN/"#""#

    if elecCount == 2:
        keyword=keyword+"2"

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

    pathStringNoStim = "dt="+str(dt)+" tStop="+str(tStop)+" p_A="+str(p_A)+" p_C="+str(p_C)+" L="+str(L)+"/"#+" myelinatedDiam="+myelDiamStr+" unmyelinatedDiam="+unmyelDiamStr
    pathString = stimulusPathString+pathStringNoStim

    suffix = {
        'elec': "electrodes/",
        'elec2': "electrodes2/",
        'draw': "draws/",
        'CAP': "CAP/",
        'CAP2': "CAP2/",
        'CAP1A': "CAPSingleAxons/",
        'CAP1A2': "CAPSingleAxons2/",
        'V': "Voltage/",
        'V2': "Voltage2/",
        'bundle': "Bundles/",
        'bundle2': "Bundles/"
    }.get(keyword,-1)

    return homeDirectory+suffix+pathString

def getFileName(recordingType, saveParams):

        directory = getDirectoryName(recordingType, **saveParams)
        if not os.path.exists(directory):
            os.makedirs(directory)

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