import numpy as np
import os
import matplotlib.pyplot as plt

locationMyelinated = '/media/carl/4ECC-1C44/Siwoo Data/Myelinated/dt=0.0025 tStop=30 pMyel=1.0 pUnmyel=0.0 L=100000 nAxons=100 (using Powley model)/bundle00000'
locationUnmyelinated = '/media/carl/4ECC-1C44/Siwoo Data/Unmyel/dt=0.0025 tStop=50 pMyel=0.0 pUnmyel=1.0 L=10000 nAxons=1000 (using Fowley model)/bundle00000'

# locationMyelinated = '/media/carl/4ECC-1C44/Siwoo Data/Myelinated/dt=0.0025 tStop=30 pMyel=1.0 pUnmyel=0.0 L=100000 nAxons=100 (using Powley model)/bundle00000'
# locationUnmyelinated = '/media/carl/4ECC-1C44/Siwoo Data/Unmyel/dt=0.0025 tStop=50 pMyel=0.0 pUnmyel=1.0 L=10000 nAxons=1000 (using Fowley model)/bundle00000'

locations = [locationUnmyelinated, locationMyelinated]

folderName = 'CAP_RecordingMechanism_recMech'
filename = 'CAP_RecordingMechanism_recMech' # .dat

startTimePlot = 1

colors = ['r', 'b']

for recMechIndex in [0,1]:
    for myelInd in [0]:
        location = locations[myelInd]

        dataLoc = os.path.join(location, folderName + str(recMechIndex), filename + str(recMechIndex) + '.dat')
        data_raw = np.transpose(np.loadtxt(dataLoc))

        # CAPraw = np.transpose(np.loadtxt(newestFile))
        time = data_raw[0, :]
        CAP = np.squeeze(data_raw[1:, :])

        plt.plot(time[time>startTimePlot], CAP[time>startTimePlot], color=colors[recMechIndex])

plt.show()