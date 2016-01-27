I assume NEURON and its depencies have been installed properly (cf NEURON website for more explanations)

****In this directory you can found:

.mod files AXNODE.mod and xtra.mod should be compiled using nrnivmodl command this create a folder which name depends if your on your OS. I attach mine (x86_64/) for Ubuntu 64 bits in case you cannot installed NEURON properly


****Files that are part of the "PyPN" package****

classes.py which contains the classes of PyPN you can instantiate among them are:
Axon, Unmyelinated, Myelinated which inherit both from Axon, Stimulus and last Bundle
Axon mainly uses copy existing from the package LFPy
Unmyelinated describes a single axon, which Hodgkin-Huxley (hh) parameters
Myelinated transcripts in Python the McIntyre myelinated axon model
Stimulus defines the stimulus signal
Bundle creates a nerve with variable proportion of myelinated and unmyelinated axons and simulate each of them separately
NOTE: In this class lies the main gain of computational efficiency, if a cluster or a machine with several CPUs and enough memory is available each axon can be run in a separate thread.


lfpcalc.py which is part of the LFPy package
refextelectrode.py which is part of the LFPy package
run_simulation.py which is part of the LFPy package

classes_plot.py which should contains the functions to save or compute the variable of interest as well as plot them nicely
(However I used "should" because while writing my dissertation I rewrote some of these functions in a specific file for each figures, and from my memories some function in classes_plot might be suboptimal such as:
for j in range(1,number_of_segments[i]):
   temp = np.loadtxt(filename, usecols=[sum(number_of_segments)-number_of_segments[i]+j+1])
   V[i][j]= temp[1:len(temp)]

which I actually rewrote, in the others scripts, using np.loadtxt only once and charging the whole file in a variable, else it is computationally inefficient.

The name of the functions should all be transparent
save_CAP_tofile, save_voltage_tofile, plot2D_CAP_fromfile, plot1D_CAP_fromfile, plotVoltage_fromfile, compute_conduction_velocity, plotConductionVelocity_fromfile, save_geometry_tofile, plot_geometry_fromfile

I however did not comment steps that should be self explanatory, i.e. most of them

****Files that are used to build a bundle and generate data*** (Equivalent of main.py)
Note: a bundle with one axon is a single axon...

bundle_CAP.py use to instantiate a Bundle in the file each variable in every instantiated dictionary has some comment associated
NOTE: the variable layout_3D which can either be DEFINE_SHAPE or PT3D refers to the corresponding NEURON function define_shape() and pt3dadd(), etc.. However, pt3d option has been added in the purpose to build a whole bundle at once without simulating each axon independently. So actually only DEFINE_SHAPE should be used.


****Files that were used to plot the figures in my dissertation****

stimulus.py plots the electrical signal shape
NB: ALL THE FOLLOWING SCRIPTS RELY ON GENERATED DATA
cv_amp_plot.py plots conduction velocity in function of the amplitude
cv_diam_plot.py plots conduction velocity in function of the diameter
outlier_plot2.py plots the outlier behavior (cf dissertation..)
plot_bipo.py plot CAP for biphasic symmetric stimulation
plot_CAP_2cm.py CAP plots for a 2cm bundle
plot_CAP2D.py plots the CAP with imshow() showing the propagation temporally and spatially
plot_CAP.py plots the CAP propagation temporally for a given position and different types of recording setting, written for monopolar stimulation
plot_CAP_bip.py same as above, some changes for biphasic symmetric stimulation
myelin_plot.py plots the AP propagation for each part of myelinated axon
unmyel_panel.py plots the AP propagation for the unmyelinated axon for different types of stimulus


****Files that were used to produce the figures in the draft****
(Equivalents of main.py)
bundle_paper.py produces data for figure1.py
bundle_paper2.py produces data for figure2.py
bundle_CAP_paper.py produces data for figure3.py

figure1.py Unmyelinated AP propagation and LFP
figure2.py Myelinated AP propagation and LFP
figure3.py spatio-temporal CAP propagation for various bundle composition
