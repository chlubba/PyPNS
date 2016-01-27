from classes_plot import *
import time

number_of_segments = 100
unmyelinated_length = 9000
diameters= [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.]
timestep = 0.0005
directory = "Voltage_data/FOR_REPORT/unmyelinated_axons_time_step0.0005axon_length9000_temperature33.0distance_stim50/"
directory = "Voltage_data/FOR_REPORT/biphasic_outlier/"
label_size = 16
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size 
text_size = 20

def loadVoltage_fromfile(filename, directory):
    filename = directory + filename
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


t0 =  time.time()
jet = plt.get_cmap('jet')
cNorm = colors.Normalize(vmin=0, vmax=number_of_segments-1)
scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)
# Using contourf to provide my colorbar info, then clearing the figure
Z = [[0,0],[0,0]]
levels = np.linspace(0,unmyelinated_length/1000.0,number_of_segments)
CS = plt.contourf(Z, levels, cmap=jet)
plt.close()
gs = gridspec.GridSpec(3, 3, width_ratios=[8,8,8], height_ratios=[6.0/16,6,6])
fig = plt.figure()
#fig,ax = plt.subplots(nrows=2, ncols=3)
amplitude = [0.2,0.3,0.7,1.2,3.4,3.5]
amplitude = [0.6,0.7,2.0,4.0,6.0,10.0]

axC = fig.add_subplot(gs[0:3])
ax = []
for k in range(2):
    for i in range(3):
        ax.append(fig.add_subplot(gs[3+k*3+i]))
        V= None
        t = None
        filename = "time_step"+str(timestep)+"unmyelinated_length9000unmyelinated_diam2.0Pulse0.01ms"+str(amplitude[i+k*3])+"nA.dat"
        filename = "p_A0.0_p_C1.0time_step0.0005recording_pos[9000]unmyelinated_length9000Pulse0.1ms"+str(amplitude[i+k*3])+"nA.dat"
        [t,V] = loadVoltage_fromfile(filename,directory)
        for j in range(number_of_segments):
            colorVal = scalarMap.to_rgba(j)
            ax[k*3+i].plot(np.array(t), np.array(V[j]), color=colorVal)
        ax[k*3+i].set_xlim([0,30])
        """if k == 0:
            ax[k*3+i].set_ylim([-70,20])"""
        #ax[k*3+i].set_ylim([-200,250])
        ax[k*3+i].set_ylim([-150,50])

        #Cvel = compute_conduction_velocity(directory+filename)
        ax[k*3+i].set_title('Stimulus amplitude '+str(amplitude[i+k*3])+' mA',x=0.52,y=0.92,fontsize =text_size+2, fontweight = 'bold',verticalalignment="bottom")

        ax[k*3+i].spines['right'].set_visible(False)
        ax[k*3+i].spines['top'].set_visible(False)
        ax[k*3+i].xaxis.set_ticks_position('none')
        ax[k*3+i].yaxis.set_ticks_position('none')
        if k == 1:
            ax[k*3+i].set_xlabel('Time [ms]', fontsize=text_size, fontweight = 'bold', labelpad=-1)
            ax[k*3+i].xaxis.set_ticks_position('bottom')
        else:
            ax[k*3+i].spines['bottom'].set_visible(False)
            for xlabel_i in ax[k*3+i].get_xticklabels():
                xlabel_i.set_visible(False)
                xlabel_i.set_fontsize(0.0)

        if i==0:
            ax[k*3+i].set_ylabel('V [mV]', fontsize=text_size, fontweight = 'bold',labelpad = -1)
            ax[k*3+i].yaxis.set_ticks_position('left')
        else:
            ax[k*3+i].spines['left'].set_visible(False)
            for ylabel_i in ax[k*3+i].get_yticklabels():
                ylabel_i.set_visible(False)
                ylabel_i.set_fontsize(0.0)
        
#divider = make_axes_locatable(ax[k*3+i])
#cax = divider.append_axes("top", size="2%", pad=0.005)
cbar =fig.colorbar(CS, cax=axC,ticks=[0, 2.5, 5, 7.5, 9],orientation="horizontal",use_gridspec=True)
cbar.set_label('Distance from electrode [mm]', rotation=0, fontsize =text_size, fontweight = 'bold', labelpad=-60)
        

#plt.tight_layout()
fig.subplots_adjust(hspace=0.1,wspace=0.1)
fig.set_size_inches(24.00,12.77)
#plt.savefig(directory+'panel.png')

print "Elapsed time for this file:" + str(time.time()-t0)
plt.show()
