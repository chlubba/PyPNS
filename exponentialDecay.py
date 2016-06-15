import numpy as np
import matplotlib.pyplot as plt

tRes = 0.1

tDelay = 10

tC = 10
aC = 5

tExp = 10
aExp = aC
cExp = -2

tD = 20
aD = -10

# def specialSignal(tRes, tDelay, tC, aC, tExp, cExp, tD, aD):
def biphasic_decaying(tRes, tDelay=5, tC=1, aC=1, tExp=1, cExp=-5, tD=2, aD=-0.2):
    """
    Defines stimulus signal shape as in <name paper>

    Args:
        tRes: time resolution of the simulation
        tDelay: length of first zero phase [ms]
        tC: length of up phase [ms]
        aC: amplitude of up phase [mA]
        tExp: length of decaying phase
        cExp: decaying coefficient
        tD: length down time
        aD: amplitude down time

    Returns:
        t: time array
        overallSignal: signal array

    """


    signal0 = np.zeros(int(tDelay / tRes))
    signal1 = np.ones(int(tC / tRes)) * aC
    tTempExp = np.arange(0, tExp, tRes)
    signal2 = np.exp(tTempExp * cExp) * aC
    signal3 = np.ones(int(tD / tRes)) * aD

    overallSignal = np.concatenate((signal0, signal1, signal2, signal3))

    t = np.arange(len(overallSignal)) * tRes

    return t, overallSignal

t, s = biphasic_decaying(tRes=0.025) # , tDelay, tC, aC, tExp, cExp, tD, aD)

plt.plot(t, s)
plt.show()