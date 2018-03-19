import numpy as np

def stationary_poisson(nsyn, lambd, tstart, tstop):
    """
    Generates nsyn stationary possion processes with rate lambda between tstart and tstop
    code from LFPy example 3

    :param nsyn: number of spike trains
    :param lambd: rate [1/s]
    :param tstart: start time [ms]
    :param tstop: stop time [ms]

    :return: nsyn spike trains stored in a matrix
    """

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


def two_phase_poisson(nsyn, lambd, tstart, tstop, cycleLength, burstiness, burstFraction=0.2):
    """
    Generates nsyn spike trains with bursting characteristic

    :param nsyn: number of spike trains
    :param lambd: rate [1/s]
    :param tstart: start time [ms]
    :param tstop: stop time [ms]
    :param cycleLength: duration on + off phase [ms]
    :param burstiness: if 0 the output is a stationary poisson process, if 1 all spikes are within the burstFraction
    :param burstFraction: share of cycleLength that where more spikes for burstiness > 0

    :return: nsyn spike trains stored in a matrix
    """

    lambd_ms = lambd * 0.001

    # number of burst-non-burst cycles
    numCycles = int(np.floor((tstop-tstart) / cycleLength) + 1)

    nonBurstLength = cycleLength * (1 - burstFraction)
    burstLength = cycleLength * burstFraction

    spiketimes = []
    for i in range(nsyn):

        # vector to save spike times into
        spikevec = np.array([])

        for j in range(numCycles):
            spikecountNonBurst = np.random.poisson(nonBurstLength * lambd_ms * (1-burstiness))
            spikecountBurst = np.random.poisson(burstLength * lambd_ms * (1 * (1-burstiness) + 1 / burstFraction * burstiness))

            spikevecNonBurst = nonBurstLength * np.random.random(spikecountNonBurst)
            spikevecBurst = nonBurstLength + burstLength * np.random.random(spikecountBurst)

            # combine both phases to one
            spikevecCycle = np.concatenate((spikevecBurst, spikevecNonBurst)) + j * cycleLength

            # append to spike train of this unit
            spikevec = np.concatenate((spikevec, spikevecCycle))

        spikevec = spikevec[spikevec < tstop] # delete spike after simulation end

        spiketimes.append(np.sort(spikevec)) #sort them too!

    return spiketimes

def generateCorrelatedSpikeTimes(nAxons, tStart=0, lambd = 1000., correlation = 0.1, tStop=300):
    """ Generate ``nAxons`` spike trains that are pairwise correlated with a factor ``correlation``. Function adapted from LFPy example 3.

    :param nAxons: number of spike streams to create
    :param tStart: start time
    :param lambd: rate [1/s]
    :param correlation: pairwise (!) correlation
    :param tStop: stop time

    :return: ``nAxons`` spike trains in a matrix

    """


    #set the numpy random seeds
    global_seed = 1234
    np.random.seed(global_seed)

    #synaptic spike times
    n_pre_syn = 1000

    #assign spike times to different units
    n_synapses = int(n_pre_syn*correlation)

    pre_syn_sptimes = stationary_poisson(nsyn=n_pre_syn, lambd=float(lambd)/n_synapses, tstart=tStart, tstop=tStop)



    signalArray = [] # np.empty(n_axons)

    for axon_id in range(nAxons):

        # re-seed the random number generator
        cell_seed = global_seed + axon_id
        np.random.seed(cell_seed)

        # Create synapse and set time of synaptic input
        pre_syn_pick = np.random.permutation(np.arange(n_pre_syn))[0:n_synapses]

        signal = []
        for i in range(n_synapses):
            signal = np.concatenate((signal, pre_syn_sptimes[pre_syn_pick[i]]))
        signal = np.sort(signal)

        signalArray.append(np.array(signal))

    return signalArray

if __name__ == "__main__":

    nAxons = 10
    spikeTrains = two_phase_poisson(nsyn=nAxons, lambd=10, tstart=0, tstop=1000, cycleLength=300,
                                                              burstiness=0.5, burstFraction=0.2)

    import matplotlib.pyplot as plt
    for i in range(nAxons):
        plt.stem(spikeTrains[i], np.ones(len(spikeTrains[i])))
    plt.show()