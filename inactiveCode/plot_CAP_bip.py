from classes_plot import *
import glob
import os
import shutil
directory = "/CAP_data2/For_report/"
path = os.getcwd() + directory 

label_size = 18
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size 
text_size = 20
amplitudes = [1.1, 1.2, 1.4, 1.6, 1.8, 2]

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

def retrieve_amp(filename):
    with open(filename, 'r') as f:
        first_line = f.readline()
    header = first_line[2:len(first_line)-1]
    bundleParameters= eval(header)
    amplitude = bundleParameters['amplitude']
    print "Amplitude: " + str(amplitude)
    return amplitude

def retrieve_time(filename):
    t = np.loadtxt(filename, unpack=True, usecols=[0])
    return t

def retrieve_mono_sum_CAP(filename):
    sum_CAP = np.loadtxt(filename, unpack=True, usecols=[1])
    return sum_CAP


def getKey(item):
    return item[0]
def getKey1(item):
    return item[1]


print path
pos = [1000+k*1000 for k in range(9)]
pos0 = [8500-1000*k for k in range(9)]
jet = plt.get_cmap('jet')
cNorm = colors.Normalize(vmin=0, vmax=len(pos)-1)
scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)
tp = []
for filename in glob.glob(os.path.join(path+"monop/",'*.dat')):
    print filename
    recording_pos = retrieve_pos(filename)
    amplitude = retrieve_amp(filename)
    tp.append([amplitude,recording_pos,filename])
    tp = sorted(tp, key=getKey)

print len(tp)
t = retrieve_time(filename)

gs = gridspec.GridSpec(2, 3, width_ratios=[8,8,8], height_ratios=[6,6])
fig = plt.figure()
ax =[]


for j in range(6):
    ax.append(fig.add_subplot(gs[j]))
    CAP = retrieve_mono_sum_CAP(tp[j][2])
    ax[j].plot(t,CAP*1000,color='b',  linestyle='-', linewidth = 2)

    ax[j].spines['right'].set_visible(False)
    ax[j].spines['top'].set_visible(False)
    ax[j].xaxis.set_ticks_position('none')
    ax[j].yaxis.set_ticks_position('none')
    ax[j].set_ylim([-150,100])
    ax[j].set_xlim([9,40])
    if j==0 or j==3:
        ax[j].yaxis.set_ticks_position('left')
        ax[j].set_ylabel('CAP ['+r'$\mu$'+'V]', fontsize=text_size, fontweight = 'bold', labelpad=0)
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
        ax[j].set_xlabel('Time [ms]', fontsize=text_size, fontweight = 'bold', labelpad=0)
    ax[j].set_title('Stimulus amplitude '+str(amplitudes[j])+' mA',x=0.52,y=0.92,fontsize =text_size+2, fontweight = 'bold',verticalalignment="bottom")
plt.suptitle('Distance pole from stimulus '+str(tp[3][1][0]/1000.0)+ ' mm / Stimulus biphasic symmetric of 0.05 ms',fontsize =text_size+3, fontweight = 'bold',verticalalignment="bottom", y = 0.94, style="italic")

fig.subplots_adjust(hspace=0.1,wspace=0.04)

plt.show()


jet = plt.get_cmap('jet')
cNorm = colors.Normalize(vmin=0, vmax=len(pos)-1)
scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)
tp = []
for filename in glob.glob(os.path.join(path+"bipo/",'*.dat')):
    print filename
    recording_pos = retrieve_pos(filename)
    amplitude = retrieve_amp(filename)
    tp.append([amplitude,recording_pos,filename])
    tp = sorted(tp, key=getKey1)
    tp = sorted(tp, key=getKey)


gs = gridspec.GridSpec(2, 3, width_ratios=[8,8,8], height_ratios=[6,6])
fig = plt.figure()
ax =[]


for j in range(6):
    CAP0 = retrieve_mono_sum_CAP(tp[j*3][2])
    CAP = retrieve_mono_sum_CAP(tp[j*3+2][2])
    CAP_bipo = CAP-CAP0

    ax.append(fig.add_subplot(gs[j]))
    ax[j].plot(t,CAP_bipo*1000,color='b',  linestyle='-', linewidth = 2)



    ax[j].spines['right'].set_visible(False)
    ax[j].spines['top'].set_visible(False)
    ax[j].xaxis.set_ticks_position('none')
    ax[j].yaxis.set_ticks_position('none')
    ax[j].set_ylim([-200,200])
    ax[j].set_xlim([9,40])
    if j==0 or j==3:
        ax[j].yaxis.set_ticks_position('left')
        ax[j].set_ylabel('CAP ['+r'$\mu$'+'V]', fontsize=text_size, fontweight = 'bold', labelpad=0)
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
        ax[j].set_xlabel('Time [ms]', fontsize=text_size, fontweight = 'bold', labelpad=0)
    ax[j].set_title('Stimulus amplitude '+str(amplitudes[j])+' mA',x=0.52,y=0.92,fontsize =text_size+2, fontweight = 'bold',verticalalignment="bottom")
plt.suptitle('First pole '+str(tp[j*3][1][0]/1000.0)+' mm from stimulus / Distance between poles '+str(tp[j*3+2][1][0]/1000.0-tp[j*3][1][0]/1000.0)+ ' mm',fontsize =text_size+3, fontweight = 'bold',verticalalignment="bottom", y = 0.94, style="italic")

fig.subplots_adjust(hspace=0.1,wspace=0.04)
plt.show()

gs = gridspec.GridSpec(2, 3, width_ratios=[8,8,8], height_ratios=[6,6])
fig = plt.figure()
ax =[]


for j in range(6):
    CAP0 = retrieve_mono_sum_CAP(tp[j*3+1][2])
    CAP = retrieve_mono_sum_CAP(tp[j*3+2][2])
    CAP_bipo = CAP-CAP0

    ax.append(fig.add_subplot(gs[j]))
    ax[j].plot(t,CAP_bipo*1000,color='b',  linestyle='-', linewidth = 2)



    ax[j].spines['right'].set_visible(False)
    ax[j].spines['top'].set_visible(False)
    ax[j].xaxis.set_ticks_position('none')
    ax[j].yaxis.set_ticks_position('none')
    ax[j].set_ylim([-200,200])
    ax[j].set_xlim([9,40])
    if j==0 or j==3:
        ax[j].yaxis.set_ticks_position('left')
        ax[j].set_ylabel('CAP ['+r'$\mu$'+'V]', fontsize=text_size, fontweight = 'bold', labelpad=0)
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
        ax[j].set_xlabel('Time [ms]', fontsize=text_size, fontweight = 'bold', labelpad=0)
    ax[j].set_title('Stimulus amplitude '+str(amplitudes[j])+' mA',x=0.52,y=0.92,fontsize =text_size+2, fontweight = 'bold',verticalalignment="bottom")
plt.suptitle('First pole '+str(tp[j*3+1][1][0]/1000.0)+' mm from stimulus / Distance between poles '+str(tp[j*3+2][1][0]/1000.0-tp[j*3+1][1][0]/1000.0)+ ' mm',fontsize =text_size+3, fontweight = 'bold',verticalalignment="bottom", y = 0.94, style="italic")

fig.subplots_adjust(hspace=0.1,wspace=0.04)
plt.show()
