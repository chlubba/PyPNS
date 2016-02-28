from classes_plot import *

import glob
import os
import shutil




label_size = 16
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size 
text_size = 20
number_of_segments = 100
unmyelinated_length = 10000

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

def loadVoltage_fromfile(filename):
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

def retrieve_parameters(filename):
    with open(filename, 'r') as f:
        first_line = f.readline()
    header = first_line[2:len(first_line)-1]
    bundleParameters= eval(header)
    myelinatedParametersA = bundleParameters['myelinated_A']
    unmyelinatedParameters = bundleParameters['unmyelinated']

    n = bundleParameters['number_elecs']
    if n != 1:
        raise TypeError("Please use plot2DCAP_fromfile(), you have recorded with more than one electrode") # handle this error in the futur with own Error class

    number_contact_points = bundleParameters['number_contact_points']
    recording_elec_pos = bundleParameters['recording_elec_pos']
    duration = bundleParameters['dur']*1e3
    return bundleParameters

def retrieve_pos(filename):
    with open(filename, 'r') as f:
        first_line = f.readline()
    header = first_line[2:len(first_line)-1]
    bundleParameters= eval(header)
    recording_elec_pos = bundleParameters['recording_elec_pos']
    print "Recording elec pos: " + str(recording_elec_pos[0])
    return recording_elec_pos

def retrieve_time(filename):
    t = np.loadtxt(filename, unpack=True, usecols=[0])
    return t

def retrieve_mono_sum_CAP(filename):
    sum_CAP = np.loadtxt(filename, unpack=True, usecols=[1])
    return sum_CAP


directory = "/FOR_PAPER/Voltage/"#"unmyelinated/"
path = os.getcwd() + directory
tp = []

for filename in glob.glob(os.path.join(path,'*.dat')):
    print filename
    nb = filename[len(path)]
    tp.append([nb,filename])
    tp = sorted(tp, key=getKey)

print "shape of tp: "
print np.shape(tp)

# pathToFile = max(glob.iglob(path+'*.[Dd][Aa][Tt]'), key=os.path.getctime)


gs = gridspec.GridSpec(3, 3, width_ratios=[8,8,8], height_ratios=[0.1,6,6])
fig = plt.figure()
axC = fig.add_subplot(gs[:3])

ax = []
for k in range(1):
    ax.append(fig.add_subplot(gs[3+k]))
    V= None
    t = None
    filename = tp[k][1]
    print filename
    [t,V] = loadVoltage_fromfile(filename)
    for j in range(number_of_segments):
        colorVal = scalarMap.to_rgba(j)
        ax[k].plot(np.array(t), np.array(V[j]), color=colorVal)
        ax[k].set_xlim([5,30])
        ax[k].set_ylim([-200,50])
    ax[k].spines['right'].set_visible(False)
    ax[k].spines['top'].set_visible(False)
    ax[k].spines['bottom'].set_visible(False)
    ax[k].xaxis.set_ticks_position('none')
    ax[k].yaxis.set_ticks_position('none')
    for xlabel_i in ax[k].get_xticklabels():
        xlabel_i.set_visible(False)
        xlabel_i.set_fontsize(0.0)
    if k==0:
        ax[k].yaxis.set_ticks_position('left')
        ax[k].set_ylabel('V [mV]', fontsize=text_size, fontweight = 'bold', labelpad=0)
    else:
        ax[k].spines['left'].set_visible(False)
        for ylabel_i in ax[k].get_yticklabels():
            ylabel_i.set_visible(False)
            ylabel_i.set_fontsize(0.0)

    
    if k == 0:
        tt = "Extracellular \n Monophasic 0.01 ms, 1 mA"
    elif k==1:
        tt = "Intracellular \n Monophasic 0.1 ms, 2 mA"
    else:
        tt = "Extracellular \n Biphasic 0.1 ms, 0.5 mA"
    ax[k].set_title(tt, y=0.82,fontsize =text_size+2, fontweight = 'bold',verticalalignment="bottom")



cbar =fig.colorbar(CS, cax=axC,ticks=[0, 2.5, 5, 7.5, 10],orientation="horizontal",use_gridspec=True)
cbar.set_label('Distance from electrode [mm]', rotation=0, fontsize =text_size, fontweight = 'bold', labelpad=-60)


# directory = "/FOR_PAPER/CAP/"#unmyelinated/"
# path = os.getcwd() + directory
# tp = []
#
# for filename in glob.glob(os.path.join(path,'*.dat')):
#     print filename
#     nb = filename[len(path)]
#     tp.append([nb,filename])
#     tp = sorted(tp, key=getKey)
#
# for k in range(3):
#     ax.append(fig.add_subplot(gs[6+k]))
#     V= None
#     t = None
#     filename = tp[k][1]
#     print filename
#     t = retrieve_time(filename)
#     V = retrieve_mono_sum_CAP(filename)
#     ax[k+3].plot(np.array(t), np.array(V*1000), color='b',  linestyle='-', linewidth = 2)
#     ax[k+3].set_xlim([5,30])
#     ax[k+3].set_ylim([-2.00,1.00])
#     ax[k+3].spines['right'].set_visible(False)
#     ax[k+3].spines['top'].set_visible(False)
#     ax[k+3].xaxis.set_ticks_position('none')
#     ax[k+3].yaxis.set_ticks_position('none')
#     if k==0:
#         ax[k+3].yaxis.set_ticks_position('left')
#         ax[k+3].set_ylabel('Extracellular potential ['+r'$\mu$'+'V]', fontsize=text_size, fontweight = 'bold', labelpad=0)
#     else:
#         ax[k+3].spines['left'].set_visible(False)
#         for ylabel_i in ax[k+3].get_yticklabels():
#             ylabel_i.set_visible(False)
#             ylabel_i.set_fontsize(0.0)
#
#
#     ax[k+3].xaxis.set_ticks_position('bottom')
#     ax[k+3].set_xlabel('Time [ms]', fontsize=text_size, fontweight = 'bold', labelpad=0)
#     if k == 0:
#         tt = "Extracellular \n Monophasic 0.01 ms, 1 mA"
#     elif k==1:
#         tt = "Intracellular \n Monophasic 0.1 ms, 2 mA"
#     else:
#         tt = "Extracellular \n Biphasic 0.1 ms, 0.5 mA"
#
#     ax[k+3].set_title(tt, y=0.82,fontsize =text_size+2, fontweight = 'bold',verticalalignment="bottom")
#

    

fig.subplots_adjust(hspace=0.1,wspace=0.08)

plt.show()
