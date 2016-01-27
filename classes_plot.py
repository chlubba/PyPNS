import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os.path
import random
import matplotlib as mpl
label_size = 20
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size 

fontP = FontProperties()
fontP.set_size(12)

def save_CAP_tofile(bundle, parameters,directory):
    filename = bundle.get_filename()
    number = 0
    filename_org = filename
    while os.path.isfile(directory+filename):
        number += 1
        print "Be careful this file name already exist ! We concatenate a number to the name to avoid erasing your previous file."
        filename = str(number) + filename_org
    print "Filename: " + filename
    filename = directory+filename
    header= repr(parameters)
    DataOut = np.array(bundle.trec)
    if bundle.number_elecs != 1:
        for i in range(len(bundle.sum_CAP)):
            DataOut = np.column_stack( (DataOut, np.array(bundle.sum_CAP[i])))
    else:        
        DataOut = np.column_stack( (DataOut, np.array(bundle.sum_CAP)))
        """
        for i in range(bundle.number_contact_points):
            DataOut = np.column_stack( (DataOut, np.array(bundle.CAP[i])))"""
    np.savetxt(filename, DataOut, header=header)

def save_voltage_tofile(bundle, parameters,directory):
    filename = bundle.get_filename()
    number = 0
    filename_org = filename
    while os.path.isfile(directory+filename):
        number += 1
        print "Be careful this file name already exist ! We concatenate a number to the name to avoid erasing your previous file."
        filename = str(number) + filename_org
    print "Filename: " + filename
    filename = directory+filename
    header= repr(parameters)
    DataOut = np.concatenate(([0],np.array(bundle.trec)))
    for w in range(bundle.number_of_axons):
        number_of_segments = 0
        try:
            number_of_segments = len(bundle.voltages[w])
        except IndexError:
            print "You have probably not set rec_v to True for this axon !"
        for i in range(number_of_segments):
            DataOut = np.column_stack( (DataOut, np.concatenate(([number_of_segments],np.array(bundle.voltages[w][i])))))
    np.savetxt(filename, DataOut, header=header)
    
def plot2D_CAP_fromfile(filename):
    with open(filename, 'r') as f:
        first_line = f.readline()
    header = first_line[2:len(first_line)-1]
    bundleParameters= eval(header)
    myelinatedParametersA = bundleParameters['myelinated_A']
    unmyelinatedParamers = bundleParameters['unmyelinated']

    n = bundleParameters['number_elecs']
    recording_elec_pos = bundleParameters['recording_elec_pos']
    duration = bundleParameters['dur']*1e3
    t = np.loadtxt(filename, unpack=True, usecols=[0])
    CAP2D = np.loadtxt(filename, unpack=True, usecols=[1])
    CAP = [ np.array(t.size) for i in range(n+1)]
    Vbrut = np.loadtxt(filename, unpack=True)
    for i in range(1,n+1):
        CAP[i-1] = Vbrut[i]
    for i in range(2,n+1):
        CAP2D = np.vstack((CAP2D,Vbrut[i]))

    X = np.linspace(0, recording_elec_pos[0], int(n))


    ### colors ###
    jet = plt.get_cmap('jet')
    cNorm = colors.Normalize(vmin=0, vmax=n-1)
    scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)
    fontP = FontProperties()
    fontP.set_size(12)

    fig =plt.figure()
    ax = plt.gca()
    CAP2D = np.clip(CAP2D*1000,-1.0*1000,0.7*1000)
    im=ax.imshow(CAP2D,cmap=cm.gist_rainbow,interpolation="none",extent=[0,duration,recording_elec_pos[0],0])
    plt.xlabel('Time('+r'$\mu$'+'s)',fontsize=24, fontweight = 'bold')
    plt.ylabel('Position along x (' +r'$\mu$'+'m)',fontsize=24, fontweight = 'bold',)
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)

    cbar =plt.colorbar(im, cax=cax)
    cbar.set_label('Amplitude CAP ('+r'$\mu$'+'V)', rotation=270, fontsize =24, fontweight = 'bold',labelpad=24)
    plt.figure()

    for i in range(n):
        colorVal = scalarMap.to_rgba(i)
        plt.plot(t,CAP[i],  color=colorVal) #label = 'CAP from ' + str(100) + ' axons at (' + str(X[i])+ ')',
    #plt.legend(loc='upper right', numpoints = 1,frameon=False , prop=fontP)
    plt.xlabel('Time(ms)')
    plt.ylabel('LFP (mV)')





def plot1D_CAP_fromfile(filename,directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
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
    t = np.loadtxt(filename, unpack=True, usecols=[0])
    sum_CAP = np.loadtxt(filename, unpack=True, usecols=[1])
    """CAP = [ np.zeros(t.size) for i in range(number_contact_points+1)]
    for i in range(1,number_contact_points+1):
        CAP[i-1] = np.loadtxt(filename, unpack=True, usecols=[i])

    radius =bundleParameters['radius_bundle']
    angles = np.linspace(0,360, number_contact_points, endpoint=False)
    Y = np.round(radius*np.cos(angles*np.pi/180.0),2)
    Z = np.round(radius*np.sin(angles*np.pi/180.0),2)

    plt.figure()
    jet = plt.get_cmap('jet')
    cNorm = colors.Normalize(vmin=0, vmax=number_contact_points)
    scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)
    for i in range(number_contact_points):
        colorVal = scalarMap.to_rgba(i)
        plt.plot(t,CAP[i]*1000, label = 'CAP from ' + str(bundleParameters['number_of_axons']) + ' axons at (' + str(Y[i]) + ',' + str(Z[i]) + ')', color=colorVal)
        plt.legend(loc='upper right', numpoints = 1,frameon=False , prop=fontP)
        plt.xlabel('Time(ms)', fontsize=24,fontweight = 'bold')
        plt.ylabel('CAP ('+r'$\mu$'+'V) at ' + str(recording_elec_pos[len(recording_elec_pos)-1])+r'$\mu$'+'m of the electrode', fontsize=24,fontweight = 'bold')"""
        

    fig =plt.figure()
    plt.plot(t,sum_CAP*1000)
    plt.xlabel('Time(ms)', fontsize=28,fontweight = 'bold')
    plt.ylabel('CAP ('+r'$\mu$'+'V) at ' + str(recording_elec_pos[len(recording_elec_pos)-1])+r'$\mu$'+'m of the electrode', fontsize=28,fontweight = 'bold')
    plt.xlim([0,bundleParameters['dur']])
    plt.ylim([-200,300])
    plt.title('Parameters: ye: '+str(bundleParameters['stim_coord'][0][1])+' duty_cycle: '+str(bundleParameters['duty_cycle'])+' amplitude: '+str(bundleParameters['amplitude']) )
    fig.set_size_inches(24.00,12.77)
    if len(recording_elec_pos) == 2:
        plt.savefig(directory+'bundle_length_'+str(unmyelinatedParameters['L'])+"recording_pos"+str(recording_elec_pos[0])+"_"+str(recording_elec_pos[1])+"p_A"+str(bundleParameters['p_A'])+"p_C"+str(bundleParameters['p_C'])+'.png')
    else:
        plt.savefig(directory+'bundle_length_'+str(unmyelinatedParameters['L'])+"recording_pos"+str(recording_elec_pos[0])+"p_A"+str(bundleParameters['p_A'])+"p_C"+str(bundleParameters['p_C'])+'amplitude_'+str(bundleParameters['amplitude'])+'.png')




def plotVoltage_fromfile(filename, directory):
    with open(filename, 'r') as f:
        first_line = f.readline()
    header = first_line[2:len(first_line)-1]
    bundleParameters= eval(header)
    myelinatedParametersA = bundleParameters['myelinated_A']
    
    unmyelinatedParameters = bundleParameters['unmyelinated']
    n = bundleParameters['number_of_axons']
    t = np.loadtxt(filename, usecols=[0])
    t = t[1:len(t)]
    #V = [[] for i in range(n)]
    V = []
    #number_of_segments = [ 0 for i in range(n)]
    number_of_segments = 0
    """for i in range(n):
        print "Currently retrieving axon number: " + str(i)
        temp = np.loadtxt(filename, usecols=[sum(number_of_segments)+1])
        number_of_segments[i] = int(temp[0])
        print "Number of segments for this axon: " + str(number_of_segments[i])
        V[i] = [np.zeros(len(t)) for k in range(number_of_segments[i])]
        V[i][0] = temp[1:len(temp)]
        for j in range(1,number_of_segments[i]):
            temp = np.loadtxt(filename, usecols=[sum(number_of_segments)-number_of_segments[i]+j+1])
            V[i][j]= temp[1:len(temp)]"""
    i = 0
    print "Currently retrieving axon number: " + str(i)
    Vbrut = np.loadtxt(filename,unpack=True)
    number_of_segments = int(Vbrut[1][0])
    print "Number of segments for this axon: " + str(number_of_segments)
    V = [np.zeros(len(t)) for k in range(number_of_segments)]
    for j in range(number_of_segments):
        temp= Vbrut[j+1]
        V[j] = temp[1:len(t)+1]
        


    Nnodes = myelinatedParametersA['Nnodes']
    nodelength = myelinatedParametersA['nodelength']
    paralength1 = myelinatedParametersA['paralength1']
    paralength2 = myelinatedParametersA['paralength2']
    interlength = myelinatedParametersA['interlength']
    for i in range(n):
        if (number_of_segments == 11*(Nnodes-1)+1): # Myelinated axon case
            # (1*nodes+2*MYSA+2*FLUT+6*STINs)*(Nnodes-1)+1node
            
            jet = plt.get_cmap('jet')
            cNorm = colors.Normalize(vmin=0, vmax=Nnodes)
            scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)
            # Using contourf to provide my colorbar info, then clearing the figure
            Z = [[0,0],[0,0]]
            levels = np.linspace(0,(nodelength*(Nnodes-1)+interlength*(Nnodes-1)*6+2*(Nnodes-1)*paralength1+ 2*(Nnodes-1)*paralength2 +0.5*nodelength)/1000.0, Nnodes-1)
            CS = plt.contourf(Z, levels, cmap=jet)
            plt.clf()
            plt.close()
            # Plot each slice as an independent subplot
            fig, ax = plt.subplots(nrows=4, ncols=1)
            for j in range(Nnodes):
                colorVal = scalarMap.to_rgba(j)
                ax[0].plot(np.array(t), np.array(V[j]), color=colorVal)
            
            ax[0].set_xlabel('Time(ms)', fontsize=20, fontweight = 'bold')
            ax[0].set_ylabel('V nodes (mV)', fontsize=20, fontweight = 'bold')
            ax[0].set_xlim([0,20])
            ax[0].set_ylim([-100,50])
            
            jet = plt.get_cmap('jet')
            cNorm = colors.Normalize(vmin=0, vmax=2*(Nnodes-1)-1)
            scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)
            for j in range(Nnodes,3*(Nnodes-1)+1):
                colorVal = scalarMap.to_rgba(j-Nnodes)
                ax[1].plot(np.array(t), np.array(V[j]), color=colorVal)

            ax[1].set_xlabel('Time(ms)', fontsize=20, fontweight = 'bold')
            ax[1].set_ylabel('V MYSA (mV)', fontsize=20, fontweight = 'bold')
            ax[1].set_xlim([0,20])
            ax[1].set_ylim([-100,50])
            
            jet = plt.get_cmap('jet')
            cNorm = colors.Normalize(vmin=0, vmax=2*(Nnodes-1)-1)
            scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)
            for j in range(3*(Nnodes-1)+1,5*(Nnodes-1)+1):
                colorVal = scalarMap.to_rgba(j-3*(Nnodes-1)+1)
                ax[2].plot(np.array(t), np.array(V[j]), color=colorVal)

            ax[2].set_xlabel('Time(ms)', fontsize=20, fontweight = 'bold')
            ax[2].set_ylabel('V FLUT (mV)', fontsize=20, fontweight = 'bold')
            ax[2].set_xlim([0,20])
            
            jet = plt.get_cmap('jet')
            cNorm = colors.Normalize(vmin=0, vmax=6*(Nnodes-1)-1)
            scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)

            for j in range(5*(Nnodes-1)+1,number_of_segments):
                colorVal = scalarMap.to_rgba(j-5*(Nnodes-1)+1)
                ax[3].plot(np.array(t), np.array(V[j]), color=colorVal)

            ax[3].set_xlabel('Time(ms)', fontsize=20, fontweight = 'bold')
            ax[3].set_ylabel('V STIN (mV)', fontsize=20, fontweight = 'bold')
            ax[3].set_xlim([0,20])
            
            plt.subplots_adjust(hspace = 0.4)
            divider = make_axes_locatable(ax[0])
            cax = divider.append_axes("top", size="10%", pad=0.4)
            #cbar =fig.colorbar(CS, cax, ticks=[0, 5, 10, 15, 20, 25, 30 ], orientation="horizontal")
            cbar =fig.colorbar(CS, cax, ticks=[0, 10,20,30,40,50, 60, 70, 80, 90, 100 ], orientation="horizontal")
            cbar.set_label('Distance from electrode (mm)', rotation=0, fontsize =19, fontweight = 'bold',labelpad=-62)
            #Cvel = compute_conduction_velocity(filename)
            #plt.title('Parameters axon fiber diameter : '+str(myelinatedParametersA['fiberD'])+', conduction velocity: '+ str(round(Cvel,2)),y=3.5, fontsize =21, fontweight = 'bold', verticalalignment="bottom")
            plt.title('Parameters axon fiber diameter : '+str(myelinatedParametersA['fiberD']),y=3.5, fontsize =21, fontweight = 'bold', verticalalignment="bottom")
            fig.set_size_inches(24.00,12.77)
            plt.savefig(directory+'axon_diam_'+str(myelinatedParametersA['fiberD'])+'.png')
            
        else:
            
            jet = plt.get_cmap('jet')
            cNorm = colors.Normalize(vmin=0, vmax=number_of_segments-1)
            scalarMap = cm.ScalarMappable(norm = cNorm, cmap= jet)
            # Using contourf to provide my colorbar info, then clearing the figure
            Z = [[0,0],[0,0]]
            levels = np.linspace(0,unmyelinatedParameters['L']/1000.0,number_of_segments)
            CS = plt.contourf(Z, levels, cmap=jet)
            plt.clf()
            plt.close()
            fig, ax = plt.subplots(nrows=1, ncols=1)
            for j in range(number_of_segments):
                colorVal = scalarMap.to_rgba(j)
                ax.plot(np.array(t), np.array(V[j]), color=colorVal)
            plt.xlabel('Time(ms)', fontsize=24, fontweight = 'bold')
            plt.ylabel('V (mV)', fontsize=24, fontweight = 'bold')
            
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("top", size="2%", pad=0.4)
            cbar =fig.colorbar(CS, cax, ticks=[0, 2.5, 5, 7.5, 9,10, 12.5, 15],orientation="horizontal")
            cbar.set_label('Distance from electrode (mm)', rotation=0, fontsize =24, fontweight = 'bold',labelpad=-62)
            #Cvel = compute_conduction_velocity(filename)
            #plt.title('Parameters axon diam :'+str(unmyelinatedParameters['diam'])+', conduction velocity: '+ str(round(Cvel,4)), y= 3.5,fontsize =22, fontweight = 'bold',verticalalignment="bottom")
            plt.title('Parameters axon diam :'+str(unmyelinatedParameters['diam']), y= 3.5,fontsize =22, fontweight = 'bold',verticalalignment="bottom")
            fig.set_size_inches(24.00,12.77)
            plt.savefig(directory+'axon_diam_'+str(unmyelinatedParameters['diam'])+'.png')
        plt.close()


def compute_conduction_velocity(filename):
    with open(filename, 'r') as f:
        first_line = f.readline()
    header = first_line[2:len(first_line)-1]
    bundleParameters= eval(header)
    myelinatedParametersA = bundleParameters['myelinated_A']
    unmyelinatedParameters = bundleParameters['unmyelinated']
    n = bundleParameters['number_of_axons'] #Actually here n should be always equal to one, because anyway only the first Cvel will be returned
    t = np.loadtxt(filename, usecols=[0])
    t = t[1:len(t)]
    V = [[] for i in range(n)]
    number_of_segments = [ 0 for i in range(n)]
    for i in range(n):
        print "Currently retrieving axon number: " + str(i)
        temp = np.loadtxt(filename, usecols=[sum(number_of_segments)+1])
        number_of_segments[i] = int(temp[0])
        print "Number of segments for this axon: " + str(number_of_segments[i])
        V[i] = [np.zeros(len(t)) for k in range(number_of_segments[i])]
        V[i][0] = temp[1:len(temp)]
        for j in range(1,number_of_segments[i]):
            temp = np.loadtxt(filename, usecols=[sum(number_of_segments)-number_of_segments[i]+j+1])
            V[i][j]= temp[1:len(temp)]

    Nnodes = myelinatedParametersA['Nnodes']
    t_max = [ [] for i in range(n)] #time to reach maximum voltage (ie AP)
    
    for i in range(n):
        if (number_of_segments[i] == 11*(Nnodes-1)+1): # Myelinated axon case
            # (1*nodes+2*MYSA+2*FLUT+6*STINs)*(Nnodes-1)+1node
            for j in range(Nnodes-1):
                t_max[i].append(t[np.argmax(V[i][j])])
                t_max[i].append(t[np.argmax(V[i][2*j+Nnodes])])
                t_max[i].append(t[np.argmax(V[i][j+2*(Nnodes-1)])])    
                for k in range(6):
                    t_max[i].append(t[np.argmax(V[i][j+Nnodes+4*(Nnodes-1)+k])])
                t_max[i].append(t[np.argmax(V[i][2*j+1+2*(Nnodes-1)])])
                t_max[i].append(t[np.argmax(V[i][2*j+1+Nnodes])])
            t_max[i].append(t[np.argmax(V[i][Nnodes])])
        else:
            for j in range(number_of_segments[i]):
                t_max[i].append(t[np.argmax(V[i][j])])


    ## Actually here either pass the type of myelinated fiber at the top of the voltage array
    # Or impose that p_A = 1
    nodelength = myelinatedParametersA['nodelength']
    paralength1 = myelinatedParametersA['paralength1']
    paralength2 = myelinatedParametersA['paralength2']
    interlength = myelinatedParametersA['interlength']
    L = unmyelinatedParameters['L']
    positions_recorded = [[] for i in range(n)]
    for k in range(n):
        if (number_of_segments[k] == 11*(Nnodes-1)+1): # myelinated axon case
            pos_nodes = []
            t_max_nodes = []
            for i in range(Nnodes):
                pos_nodes.append(nodelength*i+interlength*6*i+2*i*paralength1+ 2*i*paralength2 +0.5*nodelength)
                t_max_nodes.append(t_max[k][i+10*i])
            Cvel = float(pos_nodes[Nnodes-4]-pos_nodes[4])/(t_max_nodes[Nnodes-4]-t_max_nodes[4])
            Cvel = Cvel * 10e-6/10e-3 # m/s
        else:
            if (number_of_segments[k]%100 == 0):
                # Take an approximate of the real position recorded, because cannot retrieved it easily because we reduced number of points recorded to only 100
                for i in range(number_of_segments[k]):
                    positions_recorded[k].append(float(i+1)*L/(number_of_segments[k]))
            elif (number_of_segments[k] < 100):
                # Here take the real position recorded corresponding to the middle of segment
                for i in range(number_of_segments[k]):
                    positions_recorded[k].append(float(i)*L*3.0/(2*number_of_segments[k]))
            else:
                print "Seems like number_of_segments over 100 but not catched by myelinated case"
                print "Number of segments :" + str(number_of_segments[k])

            Cvel = float(positions_recorded[k][number_of_segments[k]-10]-positions_recorded[k][9])/(t_max[k][number_of_segments[k]-10]-t_max[k][9])
            Cvel = Cvel * 1e-6/1e-3 # m/s

        return Cvel


        
def plotConductionVelocity_fromfile(filename):
    with open(filename, 'r') as f:
        first_line = f.readline()
    header = first_line[2:len(first_line)-1]
    bundleParameters= eval(header)
    myelinatedParametersA = bundleParameters['myelinated_A']
    unmyelinatedParameters = bundleParameters['unmyelinated']
    n = bundleParameters['number_of_axons']
    t = np.loadtxt(filename, usecols=[0])
    t = t[1:len(t)]
    V = [[] for i in range(n)]
    number_of_segments = [ 0 for i in range(n)]
    for i in range(n):
        print "Currently retrieving axon number: " + str(i)
        temp = np.loadtxt(filename, usecols=[sum(number_of_segments)+1])
        number_of_segments[i] = int(temp[0])
        print "Number of segments for this axon: " + str(number_of_segments[i])
        V[i] = [np.zeros(len(t)) for k in range(number_of_segments[i])]
        V[i][0] = temp[1:len(temp)]
        for j in range(1,number_of_segments[i]):
            temp = np.loadtxt(filename, usecols=[sum(number_of_segments)-number_of_segments[i]+j+1])
            V[i][j]= temp[1:len(temp)]

    Nnodes = myelinatedParametersA['Nnodes']
    t_max = [ [] for i in range(n)] #time to reach maximum voltage (ie AP)
    v_max = [ [] for i in range(n)] # value of this maximum voltage 
    for i in range(n):
        if (number_of_segments[i] == 11*(Nnodes-1)+1): # Myelinated axon case
            # (1*nodes+2*MYSA+2*FLUT+6*STINs)*(Nnodes-1)+1node
            for j in range(Nnodes-1):
                t_max[i].append(t[np.argmax(V[i][j])])
                v_max[i].append(np.amax(V[i][j]))
                t_max[i].append(t[np.argmax(V[i][2*j+Nnodes])])
                v_max[i].append(np.amax(V[i][2*j+Nnodes]))
                t_max[i].append(t[np.argmax(V[i][j+2*(Nnodes-1)])])
                v_max[i].append(np.amax(V[i][j+2*(Nnodes-1)]))
                for k in range(6):
                    t_max[i].append(t[np.argmax(V[i][j+Nnodes+4*(Nnodes-1)+k])])
                    v_max[i].append(np.amax(V[i][j+Nnodes+4*(Nnodes-1)+k]))
                t_max[i].append(t[np.argmax(V[i][2*j+1+2*(Nnodes-1)])])
                v_max[i].append(np.amax(V[i][2*j+1+2*(Nnodes-1)]))
                t_max[i].append(t[np.argmax(V[i][2*j+1+Nnodes])])
                v_max[i].append(np.amax(V[i][2*j+1+Nnodes]))
            t_max[i].append(t[np.argmax(V[i][Nnodes])])
            v_max[i].append(np.amax(V[i][Nnodes]))
        else:
            for j in range(number_of_segments[i]):
                t_max[i].append(t[np.argmax(V[i][j])])
                v_max[i].append(np.amax(V[i][j]))


    ## Actually here either pass the type of myelinated fiber at the top of the voltage array
    # Or impose that p_A = 1
    nodelength = myelinatedParametersA['nodelength']
    paralength1 = myelinatedParametersA['paralength1']
    paralength2 = myelinatedParametersA['paralength2']
    interlength = myelinatedParametersA['interlength']
    L = unmyelinatedParameters['L']
    positions_recorded = [[] for i in range(n)]
    for k in range(n):
        if (number_of_segments[k] == 11*(Nnodes-1)+1): # myelinated axon case
            for i in range(Nnodes-1):
                positions_recorded[k].append(nodelength*i+interlength*6*i+2*i*paralength1+ 2*i*paralength2 +0.5*nodelength)
                positions_recorded[k].append(nodelength*(i+1)+interlength*6*i+2*i*paralength1+  2*i*paralength2 +0.5*paralength1)
                positions_recorded[k].append(nodelength*(i+1)+interlength*6*i+2*i*paralength1+  2*i*paralength2 + paralength1 + 0.5*paralength2)
                for j in range(6):
                    positions_recorded[k].append(nodelength*(i+1)+interlength*6*i+2*(i+1.0/2)*paralength1+  2*(i+1.0/2)*paralength2 + (0.5+j)*interlength)
                positions_recorded[k].append(nodelength*(i+1)+interlength*6*(i+1)+2*(i+1.0/2)*paralength1+  2*(i+1.0/2)*paralength2 +0.5*paralength1)
                positions_recorded[k].append(nodelength*(i+1)+interlength*6*(i+1)+2*(i+1.0/2)*paralength1+  2*(i+1.0/2)*paralength2 + paralength1 + 0.5*paralength2)
            positions_recorded[k].append(nodelength*(Nnodes-1)+interlength*(Nnodes-1)*6+2*(Nnodes-1)*paralength1+ 2*(Nnodes-1)*paralength2 +0.5*nodelength)
        else:
            if (number_of_segments[k]%100 == 0):
                # Take an approximate of the real position recorded, because cannot retrieved it easily because we reduced number of points recorded to only 100
                for i in range(number_of_segments[k]):
                    positions_recorded[k].append(float(i+1)*L/(number_of_segments[k]))
            elif (number_of_segments[k] < 100):
                # Here take the real position recorded corresponding to the middle of segment
                for i in range(number_of_segments[k]):
                    positions_recorded[k].append(float(i)*L*3.0/(2*number_of_segments[k]))
            else:
                print "Seems like number_of_segments over 100 but not catched by myelinated case"
                print "Number of segments :" + str(number_of_segments[k])

    for k in range(n):
        #print "positions:" + str(len(positions_recorded[k]))
        #print "t_max: " +str(len(t_max[k]))
        plt.figure()
        plt.title('Parameters: ye: '+str(bundleParameters['stim_coord'][0][1])+' duty_cycle: '+str(bundleParameters['duty_cycle'])+' amplitude: '+str(bundleParameters['amplitude']) )
        plt.subplot(3,1,1)
        plt.plot(positions_recorded[k], t_max[k])
        plt.xlabel('Position ('+r'$\mu$'+'m)', fontsize=20, fontweight = 'bold')
        if (number_of_segments[k] == 11*(Nnodes-1)+1): # Myelinated axon case
            plt.xlim([0,nodelength*(Nnodes-1)+interlength*(Nnodes-1)*6+2*(Nnodes-1)*paralength1+ 2*(Nnodes-1)*paralength2 +0.5*nodelength])
        plt.ylabel('Latency (ms)', fontsize=20, fontweight = 'bold')
        plt.subplot(3,1,2)
        plt.plot(positions_recorded[k], v_max[k])
        plt.xlabel('Position ('+r'$\mu$'+'m)', fontsize=20, fontweight = 'bold')
        if (number_of_segments[k] == 11*(Nnodes-1)+1): # Myelinated axon
            plt.xlim([0,nodelength*(Nnodes-1)+interlength*(Nnodes-1)*6+2*(Nnodes-1)*paralength1+ 2*(Nnodes-1)*paralength2 +0.5*nodelength])                
        plt.ylabel('V max (mV)', fontsize=20, fontweight = 'bold')
        plt.subplot(3,1,3)
        if (number_of_segments[k] == 11*(Nnodes-1)+1): # Myelinated axon case
            Cvel = np.zeros(Nnodes)
            pos_nodes = []
            t_max_nodes = []
            for i in range(Nnodes):
                pos_nodes.append(nodelength*i+interlength*6*i+2*i*paralength1+ 2*i*paralength2 +0.5*nodelength)
                t_max_nodes.append(t_max[k][i+10*i])
            for i in range(Nnodes-1):
                Cvel[i] = float(pos_nodes[i+1]-pos_nodes[i])/(t_max_nodes[i+1]-t_max_nodes[i])
            Cvel[Nnodes-1]= Cvel[Nnodes-2]
            Cvel = Cvel * 10e-6/10e-3 # m/s
            """median = np.median(Cvel)
            mean =  sum(Cvel)/Nnodes
            std = ((sum((x-median)**2 for x in Cvel))/Nnodes)**0.5
            print "Mean:" +str(mean)
            print "Std:" + str(std)
            for i in range(Nnodes):
                print abs(Cvel[i]-median)
                if abs(Cvel[i]-median)> 0.5*median:
                    Cvel[i] = median"""
            plt.plot(pos_nodes,Cvel)
            plt.xlabel('Position ('+r'$\mu$'+'m)', fontsize=20, fontweight = 'bold')
            plt.xlim([0,nodelength*(Nnodes-1)+interlength*(Nnodes-1)*6+2*(Nnodes-1)*paralength1+ 2*(Nnodes-1)*paralength2 +0.5*nodelength])                 
            plt.ylabel('Conduction Velocity (m/s)', fontsize=20, fontweight = 'bold')
            plt.ylim([0,100])
        else:
            Cvel = np.zeros(number_of_segments[k])
            for i in range(number_of_segments[k]-1):
                Cvel[i] = float(positions_recorded[k][i+1]-positions_recorded[k][i])/(t_max[k][i+1]-t_max[k][i])
            Cvel[number_of_segments[k]-1] = Cvel[number_of_segments[k]-2]
            Cvel = Cvel * 1e-6/1e-3 # m/s
            plt.plot(positions_recorded[k],Cvel)
            plt.xlabel('Position ('+r'$\mu$'+'m)', fontsize=20, fontweight = 'bold')
            plt.ylabel('Conduction Velocity (m/s)', fontsize=20, fontweight = 'bold')
            plt.ylim([0,1.5])
            
        plt.subplots_adjust(hspace = 0.3)
        print Cvel



        
def save_geometry_tofile(bundle, parameters):
    filename = "geometry"+bundle.get_filename()
    number = 0
    filename_org = filename
    while os.path.isfile("Geometry_data/"+filename):
        number += 1
        print "Be careful this file name already exist ! We concatenate a number to the name to avoid erasing your previous file."
        filename = str(number) + filename_org
    print "Filename: " + filename
    filename = "Geometry_data/"+filename
    header= repr(parameters)
    DataOut = np.array(bundle.geometry_parameters[0])
    #12 being the numbers of parameters saved by the collect_geometry function
    for i in range(1,12):
        DataOut = np.column_stack((DataOut,np.array(bundle.geometry_parameters[i])))
    np.savetxt(filename, DataOut, header=header)
    
def plot_geometry_fromfile(filename):
    with open(filename, 'r') as f:
        first_line = f.readline()
    header = first_line[2:len(first_line)-1]
    bundleParameters= eval(header)
    myelinatedParametersA = bundleParameters['myelinated_A']
    unmyelinatedParamers = bundleParameters['unmyelinated']

    n = bundleParameters['number_of_axons']
    if n != 1:
        raise TypeError("Geometry plot has been implemented only for one single axon for the moment") # handle this error in the futur with own Error class
    xstart = np.loadtxt(filename, unpack=True, usecols=[0])
    ystart = np.loadtxt(filename, unpack=True, usecols=[1])
    zstart = np.loadtxt(filename, unpack=True, usecols=[2])
    xend = np.loadtxt(filename, unpack=True, usecols=[3])
    yend = np.loadtxt(filename, unpack=True, usecols=[4])
    zend = np.loadtxt(filename, unpack=True, usecols=[5])
    area = np.loadtxt(filename, unpack=True, usecols=[6])
    diam = np.loadtxt(filename, unpack=True, usecols=[7])
    length = np.loadtxt(filename, unpack=True, usecols=[8])
    xmid = np.loadtxt(filename, unpack=True, usecols=[9])
    ymid = np.loadtxt(filename, unpack=True, usecols=[10])
    zmid = np.loadtxt(filename, unpack=True, usecols=[11])
    print xstart
    print ystart
    print zstart
    """print xmid
    print ymid
    print zmid
    print length
"""
    '''Based on LFPy example_suppl.py '''
    #creating array of points and corresponding diameters along structure
    for i in range(xend.size):
        if i == 0:
            xcoords = np.array([xmid[i]])
            ycoords = np.array([ymid[i]])
            zcoords = np.array([zmid[i]])
            diams = np.array([diam[i]])    
        else:
            if zmid[i] < 100 and zmid[i] > -100 and \
                    ymid[i] < 100 and ymid[i] > -100:
                xcoords = np.r_[xcoords, np.linspace(xstart[i],
                                            xend[i], length[i]*3)]   
                ycoords = np.r_[ycoords, np.linspace(ystart[i],
                                            yend[i], length[i]*3)]   
                zcoords = np.r_[zcoords, np.linspace(zstart[i],
                                            zend[i], length[i]*3)]   
                diams = np.r_[diams, np.linspace(diam[i], diam[i],
                                            length[i]*3)]
     #sort along depth-axis
    argsort = np.argsort(ycoords)
    """
    print xcoords
    print ycoords
    print zcoords
    print diams"""
    #plotting
    fig = plt.figure(figsize=[15, 10])
    ax = fig.add_axes([0.1, 0.1, 0.533334, 0.8], frameon=False)
    ax.scatter(xcoords[argsort], zcoords[argsort], s=diams[argsort]**2,
               c=ycoords[argsort], edgecolors='none', cmap='gray')
    #ax.plot(electrode.x, electrode.z, '.', marker='o', markersize=5, color='k')
