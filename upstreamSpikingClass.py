from neuron import h

import spikeTrainGeneration

class UpstreamSpiking(object):

    def __init__(self, nAxons, tStart, tStop, lambd = 1000., correlation = 0.1):
        self.spikeTrains = spikeTrainGeneration.generateCorrelaSpikeTimes(nAxons, tStart, tStop, lambd, correlation)
        self.axonIndex = 0 # select one specific of the nAxon spike trains in the connectAxon function

    def connectAxon(self, axon):
        # configure ExpSyn synapse
        axon.synapse = h.ExpSyn(1e-3, axon.allseclist)#1e-3
        axon.synapse.e = 10
        axon.synapse.i = 0.2
        axon.synapse.tau = 0.1

        # get spike train
        axon.spikeTrain = self.spikeTrains[self.axonIndex]
        self.axonIndex += 1

        # configure input to synapse
        axon.vecStim = h.VecStim()
        axon.spikeVec = h.Vector(axon.spikeTrain)
        # axon.spikeVec = h.Vector([1,2,3])
        axon.vecStim.play(axon.spikeVec)

        # connect synapse and VecStim input
        axon.netCon = h.NetCon(axon.vecStim, axon.synapse)
        axon.netCon.weight[0] = 1