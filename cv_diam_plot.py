from classes_plot import *
import glob
import os
label_size = 17
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size 
text_size = 20

def getKey(item):
    return item[0]



def plot_cv_diam(filename,colorVal,ax):
    diams = np.loadtxt(filename, unpack=True, usecols=[2])
    Cvels = np.loadtxt(filename, unpack=True, usecols=[3])
    tp = []
    for i in range(len(diams)):
        tp.append([diams[i],Cvels[i]])
    tp = sorted(tp, key=getKey)
    tps = [(diam, Cvel) for diam, Cvel in tp if (Cvel >= 0 and Cvel != float('Inf') and Cvel != float('-Inf'))]
    diams = []
    Cvels = []
    for i in range(len(tps)):
        diams.append(tps[i][0])
        Cvels.append(tps[i][1])        
    ax.plot(diams,Cvels,color=colorVal,  linestyle='--', marker='o', markersize = 5)


pulses = [0.001,0.002,0.004,0.01]
#pulses = [0.001,0.005,0.01,0.1]
"""gs = gridspec.GridSpec(5, 1, width_ratios=[8], height_ratios=[0.25,6,6,6,6])
fig = plt.figure()
# Plot each slice as an independent subplot
axC = fig.add_subplot(gs[0])
ax = []"""
fig, ax = plt.subplots(nrows = 2, ncols=2)
for j in range(len(pulses)):
    #directory = "/Voltage_data/unmyelinated_axons_time_step0.0005axon_length9000_temperature33.0distance_stim50/"+'Pulse'+str((pulses[j])*1.0/0.1)+'ms/'
    directory = "/Voltage_data/myelinated_axons_time_step0.0025axon_nodes21_temperature37.0stim_typeEXTRA/"+'Pulse'+str((pulses[j])*1.0/0.1)+'ms/'
    #directory = "/Voltage_data/unmyelinated_axons_time_step0.0005axon_length9000_temperature33.0stim_typeINTRA/"+'Pulse'+str((pulses[j])*1.0/0.1)+'ms/'
    folders = [ f for f in os.listdir( os.getcwd()+ directory ) if not os.path.isfile(f) ]
    folders = [f for f in folders if '#' not in f]
    tp =[]
    for i in range(len(folders)):
        tp.append([folders[i][9:len(folders[i])],folders[i]])

    tp = sorted(tp, key=getKey)
    print tp
    folders_sort = []
    for i in range(len(tp)):
        folders_sort.append(tp[i][1])
    plt.figure()
    jet = plt.get_cmap('jet')
    
    cNorm = colors.Normalize(vmin=0, vmax=len(folders))
    scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)
    # Using contourf to provide my colorbar info, then clearing the figure
    Z = [[0,0],[0,0]]
    levels = []
    for i in range(len(tp)):
        levels.append(tp[i][0])
    levels = [float(lvl) for lvl in levels]
    levels = sorted(levels)
    CS = plt.contourf(Z, levels, cmap=jet)
    plt.close()
    print levels
    for k in range(len(folders_sort)):
        path = os.getcwd() + directory +folders_sort[k] +'/'
        #print path
        for filename in glob.glob(os.path.join(path,'Cvel*.dat')):
            #print filename
            colorVal = scalarMap.to_rgba(k)
            plot_cv_diam(filename,colorVal,ax[j/2][j%2])

    if j/2 == 0:
        ax[j/2][j%2].set_ylim([0,100])
    else:
        ax[j/2][j%2].set_ylim([0,120])
    #ax[j/2][j%2].set_xlim([0,2.01])
    ax[j/2][j%2].set_xlim([5.69,16.01])
    ax[j/2][j%2].set_title('Stimulus of constant pulse length = '+str((pulses[j])*1.0/0.1)+' ms',fontsize =text_size+2, fontweight = 'bold',verticalalignment="bottom")
    
    divider = make_axes_locatable(ax[j/2][j%2])
    cax = divider.append_axes("right", size="2%", pad=0.1)

    cbar =plt.colorbar(CS, cax=cax)
    cbar.set_label('Amplitude stimulus [mA]', rotation=270, fontsize =text_size, fontweight = 'bold',labelpad=24)
    #cbar.set_label('Amplitude stimulus [nA]', rotation=270, fontsize =14, fontweight = 'bold',labelpad=24)
    
    ax[j/2][j%2].spines['right'].set_visible(False)
    ax[j/2][j%2].spines['top'].set_visible(False)
    ax[j/2][j%2].xaxis.set_ticks_position('none')
    ax[j/2][j%2].yaxis.set_ticks_position('none')
    if j==0 or j==2:
        ax[j/2][j%2].yaxis.set_ticks_position('left')
        ax[j/2][j%2].set_ylabel('Conduction velocity [m/s]', fontsize=text_size, fontweight = 'bold',labelpad = 0)
    else:
        ax[j/2][j%2].spines['left'].set_visible(False)
        for ylabel_i in ax[j/2][j%2].get_yticklabels():
            ylabel_i.set_visible(False)
            ylabel_i.set_fontsize(0.0)
        
    if j < 2:
        ax[j/2][j%2].spines['bottom'].set_visible(False)
        for xlabel_i in ax[j/2][j%2].get_xticklabels():
            xlabel_i.set_visible(False)
            xlabel_i.set_fontsize(0.0)

    else:
        ax[j/2][j%2].xaxis.set_ticks_position('bottom')
        #ax[j/2][j%2].set_xlabel('Axon diameter ['+r'$\mu$'+'m]', fontsize=text_size, fontweight = 'bold', labelpad= 0)
        ax[j/2][j%2].set_xlabel('Fibre diameter ['+r'$\mu$'+'m]', fontsize=text_size, fontweight = 'bold', labelpad= 0)
   
    
#fig.subplots_adjust(hspace=0.1,wspace=0.1)
fig.subplots_adjust(hspace=0.1,wspace=0.1)
plt.show()
