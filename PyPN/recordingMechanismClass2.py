import numpy as np
import matplotlib.pyplot as plt

class RecordingMechanism(object):

    def __init__(self, electrodePositions, extracellularPotentialMechanism, polarities=()):

        self.electrodePositions = electrodePositions

        self.numberOfPoints = np.shape(electrodePositions)[1]
        self.numberOfPoles = np.shape(electrodePositions)[2]

        if polarities == ():
            self.polarities = np.power((-1),range(self.numberOfPoles))
        else:
            self.polarities = polarities
        assert (self.numberOfPoles == len(self.polarities))

        self.extPotMech = extracellularPotentialMechanism

        self.CAP = 0
        self.CAP_axonwise = []
        # self.savePath = "" # where are the recordings stored?


    # @abstractmethod
    # def setup_recording_elec(self, bundleGuide, bundleLength):
    #     pass
    #
    # def load_one_axon(self, axonIndex):
    #
    #     directory = self.savePath
    #     filename = "electrode_"+str(axonIndex)+".dat"
    #
    #     electrodeData = np.loadtxt(os.path.join(directory, filename), unpack=True)
    #
    #     return electrodeData


    def compute_overall_CAP(self):

        arrayd_axonwise_CAPs = np.array(self.CAP_axonwise)
        self.CAP = np.sum(arrayd_axonwise_CAPs, axis=0)


    def compute_single_axon_CAP(self, axon):
        """
        1. Calculate extracellular potential (LFP) for every electrode point defined in self.electrodePositions
        2. Sum over points of single electrodes
        Args:
            axon:

        Returns: none

        """

        # The contactpoints that constitute one electrode contact (e.g. a cuff electrode ring) have to be recovered,
        # summed up together per recording location along the axon

        CAP_axonwise = np.zeros(np.shape(axon.imem)[1])

        LFPmeans = []
        # colors = ['g', 'b']
        for poleIndex in range(self.numberOfPoles):

            # calculate LFP with LFPy from membrane currents for first pole
            # LFP = self.extPotMech.calculate_LFP(axon, self.electrodePositions[:, :, poleIndex])
            LFP = self.extPotMech.calculate_LFP(axon, self.electrodePositions[:, :, poleIndex])
            CAP_axonwise += np.mean(LFP, 0)*self.polarities[poleIndex]

        #     plt.plot(np.transpose(LFP), color=colors[poleIndex])
        #
        # plt.show()

        self.CAP_axonwise.append(CAP_axonwise)



