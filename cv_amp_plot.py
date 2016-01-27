from classes_plot import *

import glob
import os
import shutil




label_size = 16
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size 
text_size = 20

diameters_m = [  5.7, 7.3, 8.7, 10.0, 11.5, 12.8, 14.0, 15.0, 16.0]
jet = plt.get_cmap('jet')
cNorm = colors.Normalize(vmin=0, vmax=20-1)#len(diameters_m)-1)#
scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)
# Using contourf to provide my colorbar info, then clearing the figure
Z = [[0,0],[0,0]]
levels = np.linspace(0.1,2.0,20)
#levels = [5.7, 7.3, 8.7, 10.0, 11.5, 12.8, 14.0, 15.0, 16.0]
CS = plt.contourf(Z, levels, cmap=jet)
plt.close()

def getKey(item):
    return item[0]
def getKey1(item):
    return item[1]


def plot_cv_amp(ax):
    nb = 0
    tps = []
    for filename in glob.glob(os.path.join(new_directory,'Cvel*.dat')):
        #print filename
        amplitude = np.loadtxt(filename, unpack=True, usecols=[1])
        diams = np.loadtxt(filename, unpack=True, usecols=[2])
        Cvels = np.loadtxt(filename, unpack=True, usecols=[3])
        tp = []
        for i in range(len(diams)):
            tp.append([round(amplitude[i],1),round(diams[i],1),round(Cvels[i],5)])
        if nb == 0:
            tps = tp
            nb += 1
        else:
            tps = np.concatenate((tps,tp), axis = 0)
        #print tps
    
    #print tps
    tps = [(round(amp,1), round(diam,1), round(Cvel,5)) for amp, diam, Cvel in tps if (Cvel >= 0 and Cvel != float('Inf') and Cvel != float('-Inf') and Cvel < 10 )]
    tps = sorted(tps, key=getKey)

    tps = sorted(tps, key=getKey1)


    for k in range(20):
    #for k in range(len(diameters_m)):
        tps_k = [(round(amp,1), round(diam,1), round(Cvel,5)) for amp, diam, Cvel in tps if (diam == round(0.1+k*0.1,1))]#round(diameters_m[k],1))]#
        #print tps_k
        colorVal = scalarMap.to_rgba(k)
        amps = []
        Cvels = []

        for i in range(len(tps_k)):
            amps.append(tps_k[i][0])
            Cvels.append(tps_k[i][2])        
        ax.plot(amps,Cvels,color=colorVal,  linestyle='--', marker='o', markersize = 5)
        
    ax.set_ylim([0,0.6])
    #ax.set_ylim([0,120])
    #EXTRA_unmyel
    ax.set_xlim([0,3.91])
    
    #INTRA
    
    #ax.set_xscale('log')
    
    #ax.set_title('Stimulus of constant pulse length = '+str((pulses[j])*1.0/0.1)+' ms',fontsize =text_size+2, fontweight = 'bold',verticalalignment="bottom")
    ax.set_title('Stimulus of constant pulse length = '+str((2*pulses[j])*1.0)+' ms',fontsize =text_size+2, fontweight = 'bold',verticalalignment="bottom")
    





def copy_cv_files_to_directory():
    
    folders = [ f for f in os.listdir( os.getcwd()+ directory ) if not os.path.isfile(f) ]
    folders = [f for f in folders if '#' not in f]
    for k in range(len(folders)):
        path = os.getcwd() + directory +folders[k] +'/'
        if not os.path.exists(new_directory):
                os.makedirs(new_directory)
        for filename in glob.glob(os.path.join(path,'Cvel*.dat')):
            shutil.copy2(filename, new_directory)

 

    
pulses = [0.001,0.002,0.004,0.01]
#pulses = [0.001,0.005,0.01,0.1]
duty = [0.005,0.01,0.025,0.05,0.025]

frequency = [0.1,0.1,0.1,0.1,0.01]
pulses = []
for j in range(len(duty)):
    pulses.append(round(duty[j]/frequency[j],2))
print pulses

gs = gridspec.GridSpec(2, 3, width_ratios=[8,8,0.25], height_ratios=[6,6])
fig = plt.figure()
axC0 = fig.add_subplot(gs[:,2])

ax = []
for j in range(len(pulses)-1):
    ax.append(fig.add_subplot(gs[j+j/2]))
    #directory = "/Voltage_data/myelinated_axons_time_step0.0025axon_nodes21_temperature37.0stim_typeEXTRA/"+'Pulse'+str((pulses[j])*1.0/0.1)+'ms/'
    #directory = "/Voltage_data/unmyelinated_axons_time_step0.0005axon_length9000_temperature33.0distance_stim50/"+'Pulse'+str((pulses[j])*1.0/0.1)+'ms/'
    directory = "/Voltage_data/Biphasic/unmyelinated_axons_time_step0.0005axon_length9000_temperature33.0distance_stim50/"+'Pulse'+str((pulses[j]))+'ms/'
    #directory = "/Voltage_data/unmyelinated_axons_time_step0.0005axon_length9000_temperature33.0stim_typeINTRA/"+'Pulse'+str((pulses[j])*1.0/0.1)+'ms/'
    new_directory = os.getcwd()+directory+"Cv_files/"
    

    copy_cv_files_to_directory()
    plot_cv_amp(ax[j])
    shutil.rmtree(new_directory)
    ax[j].spines['right'].set_visible(False)
    ax[j].spines['top'].set_visible(False)
    ax[j].xaxis.set_ticks_position('none')
    ax[j].yaxis.set_ticks_position('none')
    if j==0 or j==2:
        ax[j].yaxis.set_ticks_position('left')
        ax[j].set_ylabel('Conduction velocity [m/s]', fontsize=text_size, fontweight = 'bold', labelpad=0)
    else:
        ax[j].spines['left'].set_visible(False)
        for ylabel_i in ax[j].get_yticklabels():
            ylabel_i.set_visible(False)
            ylabel_i.set_fontsize(0.0)
        
    if j < 2:
        ax[j].spines['bottom'].set_visible(False)
        for xlabel_i in ax[j].get_xticklabels():
            xlabel_i.set_visible(False)
            xlabel_i.set_fontsize(0.0)
        #ax[j].set_ylim([60,100])

    else:
        ax[j].xaxis.set_ticks_position('bottom')
        ax[j].set_xlabel('Amplitude stimulus [mA]', fontsize=text_size, fontweight = 'bold', labelpad=0)
        #ax[j].set_xlabel('Amplitude stimulus [nA]', fontsize=text_size, fontweight = 'bold', labelpad=0)
        #ax[j].set_ylim([55,120])



cbar =fig.colorbar(CS, cax=axC0,ticks=levels,use_gridspec=True)
cbar.set_label('Axon diameter ['+r'$\mu$'+'m]', rotation=270, fontsize =text_size, fontweight = 'bold',labelpad=24)
#cbar.set_label('Fibre diameter ['+r'$\mu$'+'m]', rotation=270, fontsize =text_size, fontweight = 'bold',labelpad=24)        

   
    

fig.subplots_adjust(hspace=0.1,wspace=0.08)

plt.show()
