import matplotlib.pyplot as plt
import numpy as np
import os
import cPickle as pickle
import Tkinter as tk
import tkFileDialog as filedialog

root = tk.Tk()
root.withdraw()
filepath = filedialog.askopenfilename()

# saveDict = pickle.load(open(os.path.join('activationExtracellular', 'unmyelinated1RunsRDC0Log.dict'), "rb" )) # unmyelinated5Runs
saveDict = pickle.load(open(os.path.join(filepath), "rb" )) # unmyelinated5Runs

RDCs = saveDict['RDCs']
amplitudes = saveDict['amplitudes']
diameters = saveDict['diameters']
activationMatrix = saveDict['activationMatrix']


# filepath = filedialog.askopenfilename()
# saveDict = pickle.load(open(os.path.join(filepath), "rb" )) # unmyelinated5Runs
# activationMatrix += saveDict['activationMatrix']

print 'amplitudes:'
print amplitudes
print 'RDCs:'
print RDCs
print 'diameters:'
print diameters

(f, axarr) = plt.subplots(len(RDCs), 1)

# activationMatrix[diameterInd, RDCInd, amplitudeInd]

for RDCInd, RDC in enumerate(RDCs):

    if isinstance(axarr, np.ndarray):
        axis = axarr[RDCInd]
    else:
        axis = axarr

    cax = axis.matshow(activationMatrix[:, RDCInd, :], vmin=0, vmax=np.max(activationMatrix))
    print 'RDC %f, sum %f' % (RDC, np.sum(activationMatrix[:, RDCInd, :]))
    axis.set_xticks(range(len(amplitudes)))
    xTicks = []
    for amplitude in amplitudes:
        xTicks.append('%4.0f' % amplitude)
    axis.set_xticklabels(xTicks)
    axis.set_yticks(range(len(diameters)))
    axis.set_yticklabels(diameters)
    # axis.set_title('RDC = ' + str(RDC))

f.colorbar(cax)
plt.show()