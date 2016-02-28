from classes_plot import *

import glob
import os
import shutil




label_size = 16
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size 
text_size = 20
number_of_segments = 100
unmyelinated_length = 9000

jet = plt.get_cmap('jet')
cNorm = colors.Normalize(vmin=0, vmax=number_of_segments-1)#len(diameters_m)-1)#
scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)
# Using contourf to provide my colorbar info, then clearing the figure
Z = [[0,0],[0,0]]
levels = np.linspace(0,unmyelinated_length/1000.0,number_of_segments)
CS = plt.contourf(Z, levels, cmap=jet)
plt.close()

def getKey(item):
    return item[0]
def getKey1(item):
    return item[1]

def loadVoltage_fromfile(filename, directory):
    
    with open(filename, 'r') as f:
        first_line = f.readline()
    header = first_line[2:len(first_line)-1]
    bundleParameters= eval(header)
    myelinatedParametersA = bundleParameters['myelinated_A']
    
    unmyelinatedParameters = bundleParameters['unmyelinated']
    n = bundleParameters['number_of_axons']
    t = np.loadtxt(filename, usecols=[0])
    t = t[1:len(t)]
    V = []
    number_of_segments = 0 
    temp = np.loadtxt(filename, usecols=[1])
    number_of_segments = int(temp[0])
    print "Number of segments for this axon: " + str(number_of_segments)
    V = [np.zeros(len(t)) for k in range(number_of_segments)]
    V[0] = temp[1:len(temp)]
    Vbrute = np.loadtxt(filename, unpack=True)
    for k in range(1,number_of_segments):
        V[k]= Vbrute[k][1:len(t)+1]
    return [t,V]


directory = "/Voltage_data/FOR_REPORT/unmyelin_propag/"
path = os.getcwd() + directory 
tp = []

for filename in glob.glob(os.path.join(path,'*.dat')):
    print filename
    nb = filename[len(path)]
    tp.append([nb,filename])
    tp = sorted(tp, key=getKey)



gs = gridspec.GridSpec(2, 3, width_ratios=[8,8,8], height_ratios=[0.1,6])
fig = plt.figure()
axC = fig.add_subplot(gs[:3])

ax = []
for k in range(3):
    ax.append(fig.add_subplot(gs[3+k]))
    V= None
    t = None
    filename = tp[k][1]
    print filename
    [t,V] = loadVoltage_fromfile(filename,directory)
    for j in range(number_of_segments):
        colorVal = scalarMap.to_rgba(j)
        ax[k].plot(np.array(t), np.array(V[j]), color=colorVal)
        ax[k].set_xlim([5,30])
        ax[k].set_ylim([-140,70])
    ax[k].spines['right'].set_visible(False)
    ax[k].spines['top'].set_visible(False)
    ax[k].xaxis.set_ticks_position('none')
    ax[k].yaxis.set_ticks_position('none')
    if k==0:
        ax[k].yaxis.set_ticks_position('left')
        ax[k].set_ylabel('V [mV]', fontsize=text_size, fontweight = 'bold', labelpad=0)
    else:
        ax[k].spines['left'].set_visible(False)
        for ylabel_i in ax[k].get_yticklabels():
            ylabel_i.set_visible(False)
            ylabel_i.set_fontsize(0.0)

    
    ax[k].xaxis.set_ticks_position('bottom')
    ax[k].set_xlabel('Time [ms]', fontsize=text_size, fontweight = 'bold', labelpad=0)
    if k == 0:
        tt = "Extracellular \n Monophasic 0.01 ms, 2 mA"
    elif k==1:
        tt = "Intracellular \n Monophasic 0.1 ms, 20 mA"
    else:
        tt = "Extracellular \n Biphasic 0.1 ms, 1 mA"
    ax[k].set_title(tt, y=0.82,fontsize =text_size+2, fontweight = 'bold',verticalalignment="bottom")



cbar =fig.colorbar(CS, cax=axC,ticks=[0, 2.5, 5, 7.5, 9],orientation="horizontal",use_gridspec=True)
cbar.set_label('Distance from electrode [mm]', rotation=0, fontsize =text_size, fontweight = 'bold', labelpad=-60)
       

   
    

fig.subplots_adjust(hspace=0.1,wspace=0.08)

plt.show()
