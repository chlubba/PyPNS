from classes_plot import *
import glob
import os
import shutil
import re

text_size = 18
label_size = 14
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size 

p_A = [0,0.175,0.25,0.5,0.75,1.0]
clip_values_m = [[-0.75,0.75],[-4.6,2.8],[-7.5,3.5],[-10.5,6.7],[-17.5,11],[-35,18]]
clip_values_u = [[-0.75,0.75],[-0.6,0.4],[-1.05,0.6],[-0.83,0.6],[-0.45,0.3],[-35,18]]
number_of_segments = 100
unmyelinated_length = 9000

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


    """fig =plt.figure()
    ax = plt.gca()"""
    i = [i for i in range(len(p_A)) if (float(bundleParameters['p_A']) == p_A[i])]
    #CAP2D = np.clip(CAP2D*1000,-1.0*1000,0.7*1000)
    CAP2D = np.clip(CAP2D,clip_values_m[i[0]][0],clip_values_m[i[0]][1])
    im=ax.imshow(CAP2D,cmap=cm.gist_rainbow, interpolation="none",extent=[0,duration*1e-3,recording_elec_pos[0]*1e-3,0],aspect='auto')
    ax.set_xlim([9.8,10.5])
    ax.set_title("A-fibres "+ str(float(bundleParameters['p_A'])*100) + "%, C-fibres "+ str(float(bundleParameters['p_C'])*100)+"%",fontsize = text_size+2, fontweight = 'bold')
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.02)

    cbar =plt.colorbar(im, cax=cax)
    #cbar.set_label('Amplitude CAP ('+r'$\mu$'+'V)', rotation=270, fontsize = text_size, fontweight = 'bold',labelpad=24)
    cbar.set_label('Amplitude CAP [mV]', rotation=270, fontsize = text_size, fontweight = 'bold',labelpad=label_size+2)




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
    
    
    ax.set_xlim([9.5,12.5])
    i = [i for i in range(len(p_A)) if (float(bundleParameters['p_A']) == p_A[i])]
    ax.set_ylim([clip_values_m[i[0]][0],clip_values_m[i[0]][1]])
    ax.set_title("A-fibres "+ str(float(bundleParameters['p_A'])*100) + "%, C-fibres "+ str(float(bundleParameters['p_C'])*100)+"%",fontsize = text_size+2, fontweight = 'bold')


    
def copy_cv_files_to_directory():
    directory = "/CAP_data/2D/time_step0.0025/"
    new_directory = "CAP_data/2D/collected/"
    if os.path.exists(new_directory):
        shutil.rmtree(new_directory)
    os.makedirs(new_directory)
    path = os.getcwd() + directory
    print path
    folders = [ f for f in os.listdir( os.getcwd()+ directory ) if not os.path.isfile(f) ]
    folders = [f for f in folders if '#' not in f]
    print folders
    for k in range(len(folders)):
        folders2 = [ f for f in os.listdir( os.getcwd()+ directory+folders[k]+'/' ) if not os.path.isfile(f) ]
        folders2 = [f for f in folders2 if ('BIPOLAR' or '#') not in f]
        print folders2
        for j in range(len(folders2)):
            folders3 = [ f for f in os.listdir( os.getcwd()+ directory+folders[k]+'/' +folders2[j]+'/') if not os.path.isfile(f) ]
            folders3 = [f for f in folders3 if '9000' in f]
            print folders3
            for i in range(len(folders3)):
                path = os.getcwd() + directory +folders[k] +'/'+folders2[j]+'/'+folders3[i]+'/'
                for filename in glob.glob(os.path.join(path,'*.dat')):
                    print filename
                    shutil.copy2(filename, new_directory)
                    

def getKey(item):
    return item[0]            



#copy_cv_files_to_directory()


directory = "/CAP_data/2D/collected/"
path = os.getcwd() + directory
tps = []
for filename in glob.glob(os.path.join(path,'*.dat')):
    m = re.search('p_A(.+?)_p_C', filename)
    if m:
        found = m.group(1)
        print found
    tps.append([found,filename])
    tps = sorted(tps, key=getKey)

gs = gridspec.GridSpec(2,3, width_ratios=[8,8,8], height_ratios=[6,6])
fig = plt.figure()
ax = []
for j in range(len(p_A)):
    ax.append(fig.add_subplot(gs[j]))
    filename = tps[j][1]
    bundleParameters = load_parameters(filename)
    CAPbrut = load_CAP(filename)
    plot2D_CAP_fromfile(filename,ax[j],bundleParameters,CAPbrut)

    ax[j].spines['right'].set_visible(False)
    ax[j].spines['top'].set_visible(False)
    ax[j].xaxis.set_ticks_position('none')
    ax[j].yaxis.set_ticks_position('none')
    if j==0 or j==3:
        ax[j].yaxis.set_ticks_position('left')
        ax[j].set_ylabel('Position along x [mm]',fontsize = text_size, fontweight = 'bold')

    else:
        ax[j].spines['left'].set_visible(False)
        for ylabel_i in ax[j].get_yticklabels():
            ylabel_i.set_visible(False)
            ylabel_i.set_fontsize(0.0)
        
    if j < 3:
        ax[j].spines['bottom'].set_visible(False)
        for xlabel_i in ax[j].get_xticklabels():
            xlabel_i.set_visible(False)
            xlabel_i.set_fontsize(0.0)

    else:
        ax[j].xaxis.set_ticks_position('bottom')
        ax[j].set_xlabel('Time[ms]',fontsize = text_size, fontweight = 'bold')


fig.subplots_adjust(hspace=0.1,wspace=0.16)
plt.show()
"""


jet = plt.get_cmap('jet')
cNorm = colors.Normalize(vmin=0, vmax=number_of_segments-1)
scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)
# Using contourf to provide my colorbar info, then clearing the figure
Z = [[0,0],[0,0]]
levels = np.linspace(0,unmyelinated_length/1000.0,number_of_segments)
CS = plt.contourf(Z, levels, cmap=jet)
plt.close()
gs = gridspec.GridSpec(3,3, width_ratios=[8,8,8], height_ratios=[6.0/16,6,6])
fig = plt.figure()
ax = []
axC = fig.add_subplot(gs[0:3])
for j in range(len(p_A)):
    ax.append(fig.add_subplot(gs[j+3]))
    filename = tps[j][1]
    bundleParameters = load_parameters(filename)
    CAPbrut = load_CAP(filename)
    plot1D_CAP_fromfile(filename,ax[j],bundleParameters,CAPbrut)

    ax[j].spines['right'].set_visible(False)
    ax[j].spines['top'].set_visible(False)
    ax[j].xaxis.set_ticks_position('none')
    ax[j].yaxis.set_ticks_position('left')
    if j==0 or j==3:
        ax[j].yaxis.set_ticks_position('left')
        ax[j].set_ylabel('LFP (mV)',fontsize = text_size, fontweight = 'bold')

    #else:
     #   ax[j].spines['left'].set_visible(False)
      #  for ylabel_i in ax[j].get_yticklabels():
       #     ylabel_i.set_visible(False)
        #    ylabel_i.set_fontsize(0.0)
        
    if j < 3:
        ax[j].spines['bottom'].set_visible(False)
        for xlabel_i in ax[j].get_xticklabels():
            xlabel_i.set_visible(False)
            xlabel_i.set_fontsize(0.0)

    else:
        ax[j].xaxis.set_ticks_position('bottom')
        ax[j].set_xlabel('Time(ms)',fontsize = text_size, fontweight = 'bold')


cbar =fig.colorbar(CS, axC, ticks=[0, 2.5, 5, 7.5, 9],orientation="horizontal",use_gridspec=True)
cbar.set_label('Distance from electrode (mm)', rotation=0, fontsize =text_size, fontweight = 'bold', labelpad=-60)

fig.subplots_adjust(hspace=0.16,wspace=0.1)
plt.show()

"""
