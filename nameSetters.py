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

    if stimType in ["EXTRA", "INTRA"]:
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
        'V': "Voltage/"
    }.get(keyword,-1)

    return homeDirectory+suffix+pathString

