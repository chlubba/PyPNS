from neuron import h
from excitationMechanismClass import *
# import spikeTrainGeneration


class UpstreamSpiking(ExcitationMechanism):

    def __init__(self, spikeTrains):
        """Each axon receives a certain pattern of spikes.

        :param spikeTrains: Matrix containing the spike streams for all axons. Current problem: Axons are generated in a random manner, no way of assigning a certain spike stream to an axon of certain properties in a bundle with different axons.
        """
        self.spikeTrains = spikeTrains # spikeTrainGeneration.generateCorrelaSpikeTimes(nAxons, tStart, tStop, lambd, correlation)
        self.axonIndex = 0 # select one specific of the nAxon spike trains in the connectAxon function

        super(UpstreamSpiking, self).__init__()

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

        # add the variables to the axon instance in order to keep them alive
        # delete after simulation
        excitationMechanismVars = [synapse, vecStim, spikeVec, spikeTrain, netCon]
        axon.append_ex_mech_vars(excitationMechanismVars)

    def delete_neuron_objects(self):
        pass
