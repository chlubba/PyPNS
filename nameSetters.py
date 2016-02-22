def getDirectoryName(keyword, elecCount=2, dt=0, p_A=0, p_C=0, L=0):
    # retrieve directory name based on the purpose defined by keyword and the bundleParameters
    # 4 possible keywords:
    # elec -> temporary electrode folder
    # draw -> distrbution folder
    # CAP -> compound action potential folder
    # V -> voltage folder

    if elecCount == 2:
        keyword=keyword+"2"


    homeDirectory="/media/carl/Elements/PyPN/"#""#
    suffix = {
        'elec': "electrodes/",
        'elec2': "electrodes2/",
        'draw': "draws/",
        'CAP': "CAP_data/2D_4/time_step"+str(dt)+"/p_A"+str(p_A)+"p_C"+str(p_C)+"/RecBIPOLAR/_length"+str(L)+"/",
        'CAP2': "CAP_data/2D_4/time_step"+str(dt)+"/p_A"+str(p_A)+"p_C"+str(p_C)+"/RecMONOPOLAR/_length"+str(L)+"/",
        'V': "Voltage/"
    }.get(keyword,-1)

    return homeDirectory+suffix

