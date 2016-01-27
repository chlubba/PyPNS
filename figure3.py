from classes_plot import *
import glob
import os
import shutil
import re

text_size = 18
label_size = 14
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size 

p_A = [0,0.1,0.175,1.0]
clip_values_m = [[-0.75,0.75],[-4.6,2.8],[-4.6,2.8],[-60,40]]
#clip_values_u = [[-0.8,0.4],[-0.7,0.3],[-0.5,0.2],[-150,70]]
clip_values_u = [[-0.6,0.3],[-0.6,0.3],[-0.5,0.2],[-150,70]]
unmyelinated_length = 10000

def load_parameters(filename):
    with open(filename, 'r') as f:
        first_line = f.readline()
    header = first_line[2:len(first_line)-1]
    bundleParameters= eval(header)
    return bundleParameters

def load_CAP(filename):
    CAPbrut = np.loadtxt(filename, unpack=True)
    return CAPbrut


def plot2D_CAP_fromfile(filename,ax,bundleParameters, CAPbrut):
    myelinatedParametersA = bundleParameters['myelinated_A']
    unmyelinatedParamers = bundleParameters['unmyelinated']
    n = bundleParameters['number_elecs']
    recording_elec_pos = bundleParameters['recording_elec_pos']
    duration = bundleParameters['dur']*1e3
    t = CAPbrut[0]
    CAP = [ np.array(t.size) for i in range(n+1)]
    CAP[0] = CAPbrut[1]
    CAP2D = CAPbrut[1]
    for i in range(1,n):
        CAP[i] = CAPbrut[i+1]
        CAP2D = np.vstack((CAP2D,CAPbrut[i+1]))

    X = np.linspace(0, recording_elec_pos[0], int(n))


    ### colors ###
    jet = plt.get_cmap('jet')
    cNorm = colors.Normalize(vmin=0, vmax=n-1)
    scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)

    i = [i for i in range(len(p_A)) if (float(bundleParameters['p_A']) == p_A[i])]
    #CAP2D = np.clip(CAP2D*1000,-1.0*1000,0.7*1000)
    CAP2D = np.clip(CAP2D,clip_values_u[i[0]][0],clip_values_u[i[0]][1])
    im=ax.imshow(CAP2D,cmap=cm.gist_rainbow, interpolation="none",extent=[0,duration*1e-3,recording_elec_pos[0]*1e-3,0],aspect='auto')
    ax.set_xlim([5,30])
    #ax.set_title("A-fibres "+ str(float(bundleParameters['p_A'])*100) + "%, C-fibres "+ str(float(bundleParameters['p_C'])*100)+"%",fontsize = text_size+2, fontweight = 'bold')
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.02)
    cbar =plt.colorbar(im, cax=cax)
    if (float(bundleParameters['p_A']) ==1.0):
        ax.set_xlim([7.2,7.7])
        #cbar.set_label('Amplitude CAP ['+r'$\mu$'+'V]', rotation=270, fontsize = text_size, fontweight = 'bold',labelpad=24)
        cbar.set_label('Amplitude CAP [mV]', rotation=270, fontsize = text_size, fontweight = 'bold',labelpad=label_size+2, x=1.01)


def plot1D_CAP_fromfile(filename,ax,bundleParameters, CAPbrut):
    myelinatedParametersA = bundleParameters['myelinated_A']
    unmyelinatedParamers = bundleParameters['unmyelinated']
    n = bundleParameters['number_elecs']
    recording_elec_pos = bundleParameters['recording_elec_pos']
    duration = bundleParameters['dur']*1e3
    t = CAPbrut[0]
    CAP = [ np.array(t.size) for i in range(n+1)]
    CAP[0] = CAPbrut[1]
    CAP2D = CAPbrut[1]
    for i in range(1,n):
        CAP[i] = CAPbrut[i+1]
        CAP2D = np.vstack((CAP2D,CAPbrut[i+1]))

    X = np.linspace(0, recording_elec_pos[0], int(n))

    ### colors ###
    jet = plt.get_cmap('jet')
    cNorm = colors.Normalize(vmin=0, vmax=n-1)
    scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)
    
    for i in range(n):
        if i%10 == 0:
            colorVal = scalarMap.to_rgba(i)
            ax.plot(t,CAP[i],  color=colorVal) 
    
    
    ax.set_xlim([5,30])
    i = [i for i in range(len(p_A)) if (float(bundleParameters['p_A']) == p_A[i])]
    #ax.set_ylim([clip_values_m[i[0]][0],clip_values_m[i[0]][1]])
    ax.set_title("A-fibres "+ str(float(bundleParameters['p_A'])*100) + "%, C-fibres "+ str(float(bundleParameters['p_C'])*100)+"%",fontsize = text_size+2, fontweight = 'bold')

    
def getKey(item):
    return item[0]            


#### CAP 2D PLOTS ######

directory = "/FOR_PAPER/CAP2D/recordings/"
path = os.getcwd() + directory
tps = []
for filename in glob.glob(os.path.join(path,'*.dat')):
    m = re.search('p_A(.+?)_p_C', filename)
    if m:
        found = m.group(1)
        print found
    tps.append([found,filename])
    tps = sorted(tps, key=getKey)

gs = gridspec.GridSpec(2,4, width_ratios=[8,8,8,8], height_ratios=[7,7])
fig = plt.figure()
ax = []
for j in range(len(p_A)):
    ax.append(fig.add_subplot(gs[j+4]))
    filename = tps[j][1]
    bundleParameters = load_parameters(filename)
    CAPbrut = load_CAP(filename)
    plot2D_CAP_fromfile(filename,ax[j],bundleParameters,CAPbrut)

    ax[j].spines['right'].set_visible(False)
    ax[j].spines['top'].set_visible(False)
    ax[j].xaxis.set_ticks_position('none')
    ax[j].yaxis.set_ticks_position('none')
    if j==0:
        ax[j].yaxis.set_ticks_position('left')
        ax[j].set_ylabel('Position along x [mm]',fontsize = text_size, fontweight = 'bold')

    else:
        ax[j].spines['left'].set_visible(False)
        for ylabel_i in ax[j].get_yticklabels():
            ylabel_i.set_visible(False)
            ylabel_i.set_fontsize(0.0)
        
    if j < 0:
        ax[j].spines['bottom'].set_visible(False)
        for xlabel_i in ax[j].get_xticklabels():
            xlabel_i.set_visible(False)
            xlabel_i.set_fontsize(0.0)

    else:
        ax[j].xaxis.set_ticks_position('bottom')
        ax[j].set_xlabel('Time[ms]',fontsize = text_size, fontweight = 'bold')





### DISK DISTRIBUTION PLOTS ###
directory = "/FOR_PAPER/CAP2D/draws/"
path = os.getcwd() + directory
tps = []
for filename in glob.glob(os.path.join(path,'*.dat')):
    m = re.search('p_A(.+?)_p_C', filename)
    if m:
        found = m.group(1)
        print found
    tps.append([found,filename])
    tps = sorted(tps, key=getKey)

    
def plot_cross_section(filename,ax):
    with open(filename, 'r') as f:
        first_line = f.readline()
    header = first_line[2:len(first_line)-1]
    parameters= eval(header)
    n = parameters['number_of_axons']
    radius = parameters['radius_bundle']* np.sqrt(np.arange(n) / float(n))

    golden_angle = np.pi * (3 - np.sqrt(5))
    theta = golden_angle * np.arange(n)
    
    points = np.zeros((n, 2))
    points[:,0] = np.cos(theta)
    points[:,1] = np.sin(theta)
    points *= radius.reshape((n, 1))
    draw = (np.loadtxt(filename, unpack=True, usecols=[0]))
    diams = (np.loadtxt(filename, unpack=True, usecols=[1]))
    color = []
    for i in range(len(draw)):
        if draw[i] == 0:
            color.append('k')
        else:
            color.append('b')
    ax.scatter(points[:,1],points[:,0],c=color,s=diams*10, marker='o', edgecolors='none')
    ax.set_xlim([-175,175])
    ax.set_ylim([-175,175])
    ax.axis('equal')
    ax.set_title("A fibres "+str(parameters['p_A']*100)+"%, C fibres "+str(parameters['p_C']*100)+"%",fontsize = text_size-1, fontweight = 'bold',y = 0.92)


for j in range(len(p_A)):
    filename = tps[j][1]
    ax.append(fig.add_subplot(gs[j]))
    j= j+4
    plot_cross_section(filename,ax[j])


    ax[j].spines['right'].set_visible(False)
    ax[j].spines['top'].set_visible(False)
    ax[j].xaxis.set_ticks_position('none')
    ax[j].yaxis.set_ticks_position('none')
    if j==4:
        ax[j].yaxis.set_ticks_position('left')
        ax[j].set_ylabel('Z ['+r'$\mu$'+'m]',fontsize = text_size, fontweight = 'bold')

    else:
        ax[j].spines['left'].set_visible(False)
        for ylabel_i in ax[j].get_yticklabels():
            ylabel_i.set_visible(False)
            ylabel_i.set_fontsize(0.0)
        
    if j < 4:
        ax[j].spines['bottom'].set_visible(False)
        for xlabel_i in ax[j].get_xticklabels():
            xlabel_i.set_visible(False)
            xlabel_i.set_fontsize(0.0)

    else:
        ax[j].xaxis.set_ticks_position('bottom')
        ax[j].set_xlabel('Y ['+r'$\mu$'+'m]',fontsize = text_size, fontweight = 'bold')


fig.subplots_adjust(hspace=0.14,wspace=0.16)

plt.show()

