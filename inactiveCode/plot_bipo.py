from classes_plot import *


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

def plot_CAP_bipo(t,sum_CAP,pos0,pos1,bundleParameters):
    directory = "CAP_data2/Biphasic_bipolar_plots0.1_2/"
    if not os.path.exists(directory):
        os.makedirs(directory)
    unmyelinatedParameters = bundleParameters['unmyelinated']
    fig  =plt.figure()
    plt.plot(t,sum_CAP*1000)
    plt.xlabel('Time(ms)', fontsize=28,fontweight = 'bold')
    plt.ylabel('CAP ('+r'$\mu$'+'V) at ' + str(pos0)+'-'+str(pos1)+r'$\mu$'+'m of the electrode', fontsize=28,fontweight = 'bold')
    plt.xlim([0,bundleParameters['dur']])
    plt.title('Parameters: ye: '+str(bundleParameters['stim_coord'][0][1])+' duty_cycle: '+str(bundleParameters['duty_cycle'])+' amplitude: '+str(bundleParameters['amplitude']) )
    plt.ylim([-200,300])
    fig.set_size_inches(24.00,12.77)
    plt.savefig(directory+'bundle_length_'+str(unmyelinatedParameters['L'])+"recording_pos"+str(pos0)+'-'+str(pos1)+"p_A"+str(bundleParameters['p_A'])+"p_C"+str(bundleParameters['p_C'])+'amplitude_'+str(bundleParameters['amplitude'])+'.png')
    plt.close()


def getKey(item):
    return item[0]

import glob
import os

directory = "/CAP_data2/Biphasic_new/time_step0.0025/"
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
        folders3 = [f for f in folders3 if '#' not in f]
        print folders3
        for v in range(len(folders3)):
            folders4 = [ f for f in os.listdir( os.getcwd()+ directory+folders[k]+'/' +folders2[j]+'/'+folders3[v]+'/') if not os.path.isfile(f) ]
            folders4 = [f for f in folders4 if '#' not in f]
            print folders4
            for i in range(len(folders4)):
                print i

                path = os.getcwd() + directory +folders[k] +'/'+folders2[j]+'/'+folders3[v]+'/'+folders4[i]+'/'
                print path
                r_el_pos = []
                sum_CAP = []
                w = 0
                for filename in glob.glob(os.path.join(path,'*.dat')):
                    print filename
                    r_el_pos.append([retrieve_pos(filename),w])
                    w += 1
                    sum_CAP.append(retrieve_mono_sum_CAP(filename))
                    r_el_pos = sorted(r_el_pos, key=getKey)
                    length = len(r_el_pos)
                    increment = 0
                    t = retrieve_time(filename)
                    bundleParameters = retrieve_parameters(filename)
                    while length > 0:
                        [pos0,index0] =r_el_pos[increment]
                        for p in range(1,length):
                            [pos1,index1] =r_el_pos[increment+p]
                            sum_CAP_bipo = sum_CAP[index1]-sum_CAP[index0]
                            plot_CAP_bipo(t,sum_CAP_bipo,pos0,pos1,bundleParameters)
                        length -= 1
                        increment += 1
    
