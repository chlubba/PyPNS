from classes_plot import *
import glob
import os
import shutil
#directory = "/CAP_data2/INTRA/time_step0.0025/p_A0.2p_C0.8/MONOPOLAR/unmyelinated_length9000/"
directory = "/CAP_data2/For_report/2cm/"
path = os.getcwd() + directory 

label_size = 18
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size 
text_size = 20

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


def getKey(item):
    return item[0]


print path
pos = [1000+k*1000 for k in range(9)]
pos0 = [8500-1000*k for k in range(9)]
jet = plt.get_cmap('jet')
cNorm = colors.Normalize(vmin=0, vmax=len(pos)-1)
scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)
tp = []
for filename in glob.glob(os.path.join(path,'*.dat')):
    print filename
    recording_pos = retrieve_pos(filename)
    tp.append([recording_pos,filename])
    tp = sorted(tp, key=getKey)

print len(tp)
t = retrieve_time(filename)

gs = gridspec.GridSpec(2, 3, width_ratios=[8,8,8], height_ratios=[6,6])
fig = plt.figure()
ax =[]


for j in range(6):
    ax.append(fig.add_subplot(gs[j]))
    CAP = retrieve_mono_sum_CAP(tp[j+1][1])
    ax[j].plot(t,CAP*1000,color='b',  linestyle='-', linewidth = 2)

    ax[j].spines['right'].set_visible(False)
    ax[j].spines['top'].set_visible(False)
    ax[j].xaxis.set_ticks_position('none')
    ax[j].yaxis.set_ticks_position('none')
    ax[j].set_ylim([-1500,1000])
    ax[j].set_xlim([14,60])
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
    ax[j].set_title('Distance from stimulus '+str(tp[j+1][0][0]/1000.0)+ ' mm',fontsize =text_size+2, fontweight = 'bold',verticalalignment="bottom")

fig.subplots_adjust(hspace=0.1,wspace=0.04)

plt.show()
gs = gridspec.GridSpec(2, 2, width_ratios=[8,8], height_ratios=[6,6])
fig = plt.figure()
ax =[]


CAP5 = retrieve_mono_sum_CAP(tp[1][1])
for j in range(4):
    
    CAP = retrieve_mono_sum_CAP(tp[j+2][1])
    CAP_bipo = CAP-CAP5
    ax.append(fig.add_subplot(gs[j]))
    ax[j].plot(t,CAP_bipo*1000,color='b',  linestyle='-', linewidth = 2)

    ax[j].spines['right'].set_visible(False)
    ax[j].spines['top'].set_visible(False)
    ax[j].xaxis.set_ticks_position('none')
    ax[j].yaxis.set_ticks_position('none')
    ax[j].set_ylim([-1500,1300])
    ax[j].set_xlim([14,60])
    if j==0 or j==2:
        ax[j].yaxis.set_ticks_position('left')
        ax[j].set_ylabel('CAP ['+r'$\mu$'+'V]', fontsize=text_size, fontweight = 'bold', labelpad=0)
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
        

    else:
        ax[j].xaxis.set_ticks_position('bottom')
        ax[j].set_xlabel('Time [ms]', fontsize=text_size, fontweight = 'bold', labelpad=0)
    ax[j].set_title('Distance between poles '+str(tp[j+2][0][0]/1000.0-tp[1][0][0]/1000.0)+ ' mm',fontsize =text_size+2, fontweight = 'bold',verticalalignment="bottom")
plt.suptitle('Distance first pole from stimulus '+str(tp[1][0][0]/1000.0)+ ' mm',fontsize =text_size+3, fontweight = 'bold',verticalalignment="bottom", y = 0.94, style="italic")
fig.subplots_adjust(hspace=0.1,wspace=0.04)

plt.show()

"""gs1 = gridspec.GridSpec(2, 2, width_ratios=[8,8], height_ratios=[6,6])
gs1.update(left=0.05, right=0.49, wspace=0.04,hspace=0.1)"""
fig = plt.figure()
ax =[]
for j in range(4):
    CAP0 = retrieve_mono_sum_CAP(tp[j+2-1][1])
    CAP = retrieve_mono_sum_CAP(tp[j+2+1][1])
    CAP_bipo = CAP-CAP0

    ax.append(fig.add_subplot(gs[j]))
    ax[j].plot(t,CAP_bipo*1000,color='b',  linestyle='-', linewidth = 2)

    ax[j].spines['right'].set_visible(False)
    ax[j].spines['top'].set_visible(False)
    ax[j].xaxis.set_ticks_position('none')
    ax[j].yaxis.set_ticks_position('none')
    ax[j].set_ylim([-1500,1300])
    ax[j].set_xlim([14,60])
    if j==0 or j==2:
        ax[j].yaxis.set_ticks_position('left')
        ax[j].set_ylabel('CAP ['+r'$\mu$'+'V]', fontsize=text_size, fontweight = 'bold', labelpad=0)
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
        

    else:
        ax[j].xaxis.set_ticks_position('bottom')
        ax[j].set_xlabel('Time [ms]', fontsize=text_size, fontweight = 'bold', labelpad=0)
    ax[j].set_title('Distance first pole from stimulus '+str(tp[j+2-1][0][0]/1000.0)+ ' mm',fontsize =text_size+2, fontweight = 'bold',verticalalignment="bottom")
plt.suptitle('Distance between poles '+str(tp[j+2+1][0][0]/1000.0-tp[j+2-1][0][0]/1000.0)+ ' mm',fontsize =text_size+3, fontweight = 'bold',verticalalignment="bottom", y = 0.94, style="italic")

fig.subplots_adjust(hspace=0.1,wspace=0.04)
fig = plt.figure()
"""gs2 = gridspec.GridSpec(2, 2, width_ratios=[8,8], height_ratios=[6,6])
gs2.update(left=0.51, right=0.95, wspace = 0.04, hspace=0.1)"""
ax =[]
for j in range(4):
    CAP0 = retrieve_mono_sum_CAP(tp[j+1-1][1])
    CAP = retrieve_mono_sum_CAP(tp[j+1+2][1])
    CAP_bipo = CAP-CAP0

    ax.append(fig.add_subplot(gs[j]))
    ax[j].plot(t,CAP_bipo*1000,color='b',  linestyle='-', linewidth = 2)

    ax[j].spines['right'].set_visible(False)
    ax[j].spines['top'].set_visible(False)
    ax[j].xaxis.set_ticks_position('none')
    ax[j].yaxis.set_ticks_position('none')
    ax[j].set_ylim([-1500,1300])
    ax[j].set_xlim([14,60])
    if j==0 or j==2:
    #if False:
        ax[j].yaxis.set_ticks_position('left')
        ax[j].set_ylabel('CAP ['+r'$\mu$'+'V]', fontsize=text_size, fontweight = 'bold', labelpad=0)
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
        

    else:
        ax[j].xaxis.set_ticks_position('bottom')
        ax[j].set_xlabel('Time [ms]', fontsize=text_size, fontweight = 'bold', labelpad=0)
    ax[j].set_title('Distance first pole from stimulus '+str(tp[j+1-1][0][0]/1000.0)+ ' mm',fontsize =text_size+2, fontweight = 'bold',verticalalignment="bottom")
plt.suptitle('Distance between poles '+str(tp[j+1+2][0][0]/1000.0-tp[j+1-1][0][0]/1000.0)+ ' mm',fontsize =text_size+3, fontweight = 'bold',verticalalignment="bottom", y = 0.94, style="italic")
fig.subplots_adjust(hspace=0.1,wspace=0.04)

plt.show()


