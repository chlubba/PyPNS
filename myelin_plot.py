from classes_plot import *


Nnodes = 21
diameters= [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.]
timestep = 0.0005
directory = "Voltage_data/FOR_REPORT/myelinated_axons_time_step"+str(timestep)+"axon_nodes"+str(Nnodes)+"_temperature"+str(37.0)+"distance_stim"+str(50)+"/"
filename = "myelinated_nodes21myelinated_diam11.5Pulse0.01ms1.0nA.dat"
label_size = 16
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size 
text_size = 18

def loadVoltage_fromfile(filename, directory):
    filename = directory + filename
    t = np.loadtxt(filename, usecols=[0])
    t = t[1:len(t)]
    V = []
    number_of_segments = 0 
    temp = np.loadtxt(filename, usecols=[1])
    number_of_segments = int(temp[0])
    print "Number of segments for this axon: " + str(number_of_segments)
    V = [np.zeros(len(t)) for k in range(number_of_segments)]
    V[0] = temp[1:len(temp)]
    Vbrut = np.loadtxt(filename, unpack=True)
    for j in range(1,number_of_segments):
        temp = Vbrut[j+1]
        V[j]= temp[1:len(temp)]
    return [t,V]

def loadparameters_fromfile(filename, directory):
    filename = directory + filename
    with open(filename, 'r') as f:
        first_line = f.readline()
    header = first_line[2:len(first_line)-1]
    bundleParameters= eval(header)
    myelinatedParametersA = bundleParameters['myelinated_A']
    myelinatedParametersB = bundleParameters['myelinated_B']
    unmyelinatedParameters = bundleParameters['unmyelinated']

    return myelinatedParametersA

def compute_cv2(t,V,myelinatedParametersA):
    t_max = [] #time to reach maximum voltage (ie AP)
    Nnodes = myelinatedParametersA['Nnodes']
    nodelength = myelinatedParametersA['nodelength']
    paralength1 = myelinatedParametersA['paralength1']
    paralength2 = myelinatedParametersA['paralength2']
    interlength = myelinatedParametersA['interlength']

    for j in range(Nnodes):
        t_max.append(t[np.argmax(V[j])])
    pos_nodes = []
    t_max_nodes = []
    for i in range(Nnodes):
        pos_nodes.append(nodelength*i+interlength*6*i+2*i*paralength1+ 2*i*paralength2 +0.5*nodelength)
        
    Cvel = float(pos_nodes[Nnodes-4]-pos_nodes[4])/(t_max[Nnodes-4]-t_max[4])
    Cvel = Cvel * 1e-6/1e-3 # m/s
    return round(Cvel,2)


myelinatedParametersA = loadparameters_fromfile(filename, directory)
Nnodes = myelinatedParametersA['Nnodes']
nodelength = myelinatedParametersA['nodelength']
paralength1 = myelinatedParametersA['paralength1']
paralength2 = myelinatedParametersA['paralength2']
interlength = myelinatedParametersA['interlength']
[t,V] = loadVoltage_fromfile(filename,directory)
Cvel = compute_cv2(t,V,myelinatedParametersA)

jet = plt.get_cmap('jet')
cNorm = colors.Normalize(vmin=0, vmax=Nnodes)
scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)
# Using contourf to provide my colorbar info, then clearing the figure
Z = [[0,0],[0,0]]
levels = np.linspace(0,(nodelength*(Nnodes-1)+interlength*(Nnodes-1)*6+2*(Nnodes-1)*paralength1+ 2*(Nnodes-1)*paralength2 +0.5*nodelength)/1000.0, Nnodes-1)
CS = plt.contourf(Z, levels, cmap=jet)
plt.clf()
plt.close()
gs = gridspec.GridSpec(5, 1, width_ratios=[8], height_ratios=[0.5,6,6,6,6])
fig = plt.figure()
# Plot each slice as an independent subplot
axC = fig.add_subplot(gs[0])
ax = []
ax.append(fig.add_subplot(gs[1]))
for j in range(Nnodes):
    colorVal = scalarMap.to_rgba(j)
    ax[0].plot(np.array(t), np.array(V[j]), color=colorVal)
ax[0].set_ylabel('V nodes [mV]', fontsize=text_size, fontweight = 'bold')
ax[0].set_ylim([-100,50])

jet = plt.get_cmap('jet')
cNorm = colors.Normalize(vmin=0, vmax=2*(Nnodes-1)-1)
scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)
ax.append(fig.add_subplot(gs[2]))
for j in range(Nnodes,3*(Nnodes-1)+1):
    colorVal = scalarMap.to_rgba(j-Nnodes)
    ax[1].plot(np.array(t), np.array(V[j]), color=colorVal)
ax[1].set_ylabel('V MYSA [mV]', fontsize=text_size, fontweight = 'bold')
ax[1].set_ylim([-100,50])

jet = plt.get_cmap('jet')
cNorm = colors.Normalize(vmin=0, vmax=2*(Nnodes-1)-1)
scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)
ax.append(fig.add_subplot(gs[3]))
for j in range(3*(Nnodes-1)+1,5*(Nnodes-1)+1):
    colorVal = scalarMap.to_rgba(j-3*(Nnodes-1)+1)
    ax[2].plot(np.array(t), np.array(V[j]), color=colorVal)
ax[2].set_ylabel('V FLUT [mV]', fontsize=text_size, fontweight = 'bold')

jet = plt.get_cmap('jet')
cNorm = colors.Normalize(vmin=0, vmax=6*(Nnodes-1)-1)
scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)
ax.append(fig.add_subplot(gs[4]))
for j in range(5*(Nnodes-1)+1,(1+2+2+6)*(Nnodes-1)+1):
    colorVal = scalarMap.to_rgba(j-5*(Nnodes-1)+1)
    ax[3].plot(np.array(t), np.array(V[j]), color=colorVal)
ax[3].set_ylabel('V STIN [mV]', fontsize=text_size, fontweight = 'bold')


cbar =fig.colorbar(CS, cax=axC,ticks=[j*5 for j in range(21)],orientation="horizontal",use_gridspec=True)
cbar.set_label('Distance from electrode [mm]', rotation=0, fontsize =text_size, fontweight = 'bold', labelpad=-60)

for j in range(4):
    ax[j].spines['right'].set_visible(False)
    ax[j].spines['top'].set_visible(False)
    ax[j].xaxis.set_ticks_position('none')
    ax[j].yaxis.set_ticks_position('left')
    ax[j].set_xlim([0,20])
    if j == 3:
        ax[j].set_xlabel('Time [ms]', fontsize=text_size, fontweight = 'bold')
        
    else:
        ax[j].spines['bottom'].set_visible(False)
        for xlabel_i in ax[j].get_xticklabels():
            xlabel_i.set_visible(False)
            xlabel_i.set_fontsize(0.0)


axC.set_title('Axon fibre diameter '+str(myelinatedParametersA['fiberD'])+r' $\mu$'+'m'+", Conduction velocity "+str(Cvel)+" m/s" ,y=3.5, fontsize =text_size+2, fontweight = 'bold', verticalalignment="bottom")

            

fig.subplots_adjust(hspace=0.15,wspace=0.1)
fig.set_size_inches(24.00,12.77)

plt.show()
