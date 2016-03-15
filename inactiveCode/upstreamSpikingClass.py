from neuron import h
from ExcitationMechanism import *

import spikeTrainGeneration

class UpstreamSpiking(ExcitationMechanism):

    def __init__(self, nAxons, tStart, tStop, lambd = 1000., correlation = 0.1):
        self.spikeTrains = spikeTrainGeneration.generateCorrelaSpikeTimes(nAxons, tStart, tStop, lambd, correlation)
        self.axonIndex = 0 # select one specific of the nAxon spike trains in the connectAxon function

    def connect_axon(self, axon):

        # configure ExpSyn synapse
        synapse = h.ExpSyn(1e-3, axon.allseclist)#1e-3
        synapse.e = 10
        synapse.i = 0.2
        synapse.tau = 0.1

        # get spike train
        spikeTrain = self.spikeTrains[self.axonIndex]
        self.axonIndex += 1

        # configure input to synapse
        vecStim = h.VecStim()
        spikeVec = h.Vector(spikeTrain)
        # axon.spikeVec = h.Vector([1,2,3])
        vecStim.play(spikeVec)

        # connect synapse and VecStim input
        netCon = h.NetCon(vecStim, synapse)
        netCon.weight[0] = 1

        excitationMechanismVars = [synapse, vecStim, spikeVec, spikeTrain, netCon]

        axon.appendExMechVars(excitationMechanismVars)
