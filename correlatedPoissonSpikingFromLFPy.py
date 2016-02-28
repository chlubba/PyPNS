#!/usr/bin/env python
'''
LFPs from a population of cells relying on MPI
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
import os
import sys
if sys.version < '3':
    from urllib2 import urlopen
else:
    from urllib.request import urlopen
import zipfile
import LFPy
import neuron



def stationary_poisson(nsyn,lambd,tstart,tstop):
    ''' Generates nsyn stationary possion processes with rate lambda between tstart and tstop'''
    # code from LFPy example 3
    interval_s = (tstop-tstart)*.001
    spiketimes = []
    for i in range(nsyn):
        spikecount = np.random.poisson(interval_s*lambd)
        spikevec = np.empty(spikecount)
        if spikecount==0:
            spiketimes.append(spikevec)
        else:
            spikevec = tstart + (tstop-tstart)*np.random.random(spikecount)
            spiketimes.append(np.sort(spikevec)) #sort them too!

    return spiketimes


wanted_correlation = 0.1

#number of units
n_axons = 2
axon_id = 1

#set the numpy random seeds
global_seed = 1234
np.random.seed(global_seed)

#synaptic spike times
n_pre_syn = 10000
pre_syn_sptimes = stationary_poisson(nsyn=n_pre_syn, lambd=5., tstart=0, tstop=300)

#assign spike times to different units
n_synapses = int(n_pre_syn*wanted_correlation)


# FIRST AXON

# re-seed the random number generator
cell_seed = global_seed + axon_id
np.random.seed(cell_seed)

# Create synapse and set time of synaptic input
pre_syn_pick = np.random.permutation(np.arange(n_pre_syn))[0:n_synapses]

signal1 = []
for i in range(n_synapses):
    signal1 = np.concatenate((signal1, pre_syn_sptimes[pre_syn_pick[i]])) #[1:n_synapses]

# SECOND AXON

axon_id = 2

# re-seed the random number generator
cell_seed = global_seed + axon_id
np.random.seed(cell_seed)

# Create synapse and set time of synaptic input
pre_syn_pick = np.random.permutation(np.arange(n_pre_syn))[0:n_synapses]

signal2 = []
for i in range(n_synapses):
    signal2 = np.concatenate((signal2, pre_syn_sptimes[pre_syn_pick[i]])) #[1:n_synapses]


sigHist1 = np.histogram(signal1, bins = range(301))[0]
sigHist2 = np.histogram(signal2, bins = range(301))[0]

print 'correlation = ' + str(np.corrcoef(sigHist1, sigHist2))


