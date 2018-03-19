import os
import sys
import glob
import cPickle as pickle

def get_bundle_directory(paramDict, new = False):
    """Create the directory where all output of the simulation is saved in. If no ``saveLocation`` is specified in the ``paramDict``, the working directory is used as a basis. For each simulation an individual subfolder will be generated.

    :param paramDict: dictionary containing all information charaterizing the ``Bundle`` and the wanted save location.
    :param new: True if new directory is to be generated. False if directory of previous calculation wants to be retrieved.

    :return: location for simulation output
    :rtype: string

    """

    # read out dictionary of parameters (more elegant method possible?)
    dt=paramDict['timeRes']
    tStop = paramDict['tStop']
    pMyel = paramDict['pMyel']
    L = paramDict['length']
    numberOfAxons = paramDict['numberOfAxons']

    paramsMyel = paramDict['paramsMyel']
    paramsUnmyel = paramDict['paramsUnmyel']
    myelinatedDiam = paramsMyel['fiberD']
    unmyelinatedDiam = paramsUnmyel['fiberD']

    homeDirectory = paramDict['saveLocation']

    pUnmyel = 1 - pMyel

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
    Get directory name based on the purpose defined by keyword and the ``basePath`` where all output will be stored.

    :param keyword: prefix of the file name
    :param basePath: bundle basepath

    :return: directory name

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
    """Get the file name for certain output data including the folder.

    :param recordingType: Type of output
    :param basePath: base path where all output is stored
    :param newFile: If true, new unambiguous file name will be returned (for saving) otherwise most recent existing filename is given (loading).
    :param directoryType: directory name can be specified. If unset, directory name is equal to file prefix.

    :return: file name includig folder

    """

    if isinstance(directoryType, bool):
        directory = get_directory_name(recordingType, basePath)
    else:
        directory = get_directory_name(directoryType, basePath)

    # filename = 'recording.dat'
    filename = recordingType+'.npy'

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
    """Use ``pickle`` to save the ``Bundle``.

    :param bundle: ``Bundle`` class
    """

    # save the bundle definition file
    bundleSaveLocation = bundle.basePath
    pickle.dump(bundle,open( os.path.join(bundleSaveLocation, 'bundle.cl'), "wb" ))

def open_recent_bundle(parameters):
    """Open most recent bundle in at the specified base path with save parameters as specified. Not all characteristics of ``Bundle`` are captured in the save parameters so beware of loading the wrong ``Bundle``.

    :param parameters: bundle parameters as a dictionary

    :return: ``Bundle`` instance
    """
    bundleSaveLocation = get_bundle_directory(parameters, new=False)
    try:
        bundle = pickle.load(open(os.path.join(bundleSaveLocation, 'bundle.cl'), "rb" ))
        return bundle
    except:
        print 'No bundle with these parameters has been generated yet. Set calculationFlag to True.'


def open_bundle_from_location(bundleSaveLocation):
    """Like ``open_recent_bundle`` but with a location as an argument. More reliable.

    :param bundleSaveLocation: Exact location of the pickled bundle file.

    :return: ``Bundle`` instance
    """

    try:
        bundle = pickle.load(open(os.path.join(bundleSaveLocation, 'bundle.cl'), "rb" ))
        return bundle
    except:
        print("Unexpected error:", sys.exc_info()[0])
        print 'No bundle with these parameters has been generated yet. Set calculationFlag to True.'
        raise

