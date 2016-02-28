from classes_plot import *
import glob

Nnodes = 41

timestep = 0.0025
directory = "/FOR_PAPER/Voltage/myelinated/"
path = os.getcwd() + directory
for filename in glob.glob(os.path.join(path,'*.dat')):
    print filename
    filename = filename

label_size = 16
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size 
text_size = 18

def getKey(item):
    return item[0]

def loadVoltage_fromfile(filename):
    
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

def loadparameters_fromfile(filename):
    
    with open(filename, 'r') as f:
        first_line = f.readline()
    header = first_line[2:len(first_line)-1]
    bundleParameters= eval(header)
    myelinatedParametersA = bundleParameters['myelinated_A']
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


myelinatedParametersA = loadparameters_fromfile(filename)
Nnodes = myelinatedParametersA['Nnodes']
nodelength = myelinatedParametersA['nodelength']
paralength1 = myelinatedParametersA['paralength1']
paralength2 = myelinatedParametersA['paralength2']
interlength = myelinatedParametersA['interlength']
[t,V] = loadVoltage_fromfile(filename)
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
gs = gridspec.GridSpec(5, 2, width_ratios=[8,8], height_ratios=[0.5,6,6,6,6])
fig = plt.figure()
# Plot each slice as an independent subplot
axC = fig.add_subplot(gs[0])
ax = []
ax.append(fig.add_subplot(gs[1*2]))
for j in range(Nnodes):
    colorVal = scalarMap.to_rgba(j)
    ax[0].plot(np.array(t), np.array(V[j]), color=colorVal)
ax[0].set_ylabel('V nodes [mV]', fontsize=text_size, fontweight = 'bold')
ax[0].set_ylim([-100,50])

jet = plt.get_cmap('jet')
cNorm = colors.Normalize(vmin=0, vmax=2*(Nnodes-1)-1)
scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)
ax.append(fig.add_subplot(gs[2*2]))
for j in range(Nnodes,3*(Nnodes-1)+1):
    colorVal = scalarMap.to_rgba(j-Nnodes)
    ax[1].plot(np.array(t), np.array(V[j]), color=colorVal)
ax[1].set_ylabel('V MYSA [mV]', fontsize=text_size, fontweight = 'bold')
ax[1].set_ylim([-100,50])

jet = plt.get_cmap('jet')
cNorm = colors.Normalize(vmin=0, vmax=2*(Nnodes-1)-1)
scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)
ax.append(fig.add_subplot(gs[2*3]))
for j in range(3*(Nnodes-1)+1,5*(Nnodes-1)+1):
    colorVal = scalarMap.to_rgba(j-3*(Nnodes-1)+1)
    ax[2].plot(np.array(t), np.array(V[j]), color=colorVal)
ax[2].set_ylabel('V FLUT [mV]', fontsize=text_size, fontweight = 'bold')

jet = plt.get_cmap('jet')
cNorm = colors.Normalize(vmin=0, vmax=6*(Nnodes-1)-1)
scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)
ax.append(fig.add_subplot(gs[2*4]))
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



directory = "/FOR_PAPER/CAP/myelinated"
path = os.getcwd() + directory 
tp = []

for filename in glob.glob(os.path.join(path,'*.dat')):
    print filename
    pos = retrieve_pos(filename)
    tp.append([pos,filename])
    tp = sorted(tp, key=getKey)
    
for k in range(4):
    ax.append(fig.add_subplot(gs[1+(k+1)*2]))
    V= None
    t = None
    filename = tp[k][1]
    print filename
    t = retrieve_time(filename)
    V = retrieve_mono_sum_CAP(filename)
    ax[k+4].plot(np.array(t), np.array(V*1000), color='b',  linestyle='-', linewidth = 2)
    ax[k+4].set_xlim([0,20])
    ax[k+4].set_ylim([-20.00,10.00])
    ax[k+4].spines['right'].set_visible(False)
    ax[k+4].spines['top'].set_visible(False)
    ax[k+4].xaxis.set_ticks_position('none')
    ax[k+4].yaxis.set_ticks_position('none')
    if True:
        ax[k+4].yaxis.set_ticks_position('left')
        if k == 0:
            section = "nodes"
        elif k == 1:
            section = "FLUT"
        elif k == 2:
            section = "MYSA"
        else:
            section = "STIN"
        ax[k+4].set_ylabel("Vext "+section+' ['+r'$\mu$'+'V]', fontsize=text_size, fontweight = 'bold', labelpad=0)
    else:
        ax[k+4].spines['left'].set_visible(False)
        for ylabel_i in ax[k+4].get_yticklabels():
            ylabel_i.set_visible(False)
            ylabel_i.set_fontsize(0.0)
    if k < 3:
        ax[k+4].spines['bottom'].set_visible(False)
        for xlabel_i in ax[k+4].get_xticklabels():
            xlabel_i.set_visible(False)
            xlabel_i.set_fontsize(0.0)

    else:
        ax[k+4].xaxis.set_ticks_position('bottom')
        ax[k+4].set_xlabel('Time(ms)',fontsize = text_size, fontweight = 'bold')


#axC.set_title('Axon fibre diameter '+str(myelinatedParametersA['fiberD'])+r' $\mu$'+'m'+", Conduction velocity "+str(Cvel)+" m/s" ,y=3.5, fontsize =text_size+2, fontweight = 'bold', verticalalignment="bottom")
            

fig.subplots_adjust(hspace=0.15,wspace=0.1)
fig.set_size_inches(24.00,12.77)
plt.show()
