import numpy as np
import matplotlib.pyplot as plt

class RecordingMechanism(object):

    def __init__(self, electrodePositions, extracellularPotentialMechanism, polarities=()):
        """``RecordingMechanism`` calculates the extracellular potental as picked up by an electrode.

        :param electrodePositions: Positions of the recording electrodes. 3D numpy array for (electrode points) x (position in space of each electrode point ) x (poles)
        :param extracellularPotentialMechanism: ``ExtracellularMechanism`` to convert membrane current into voltage
        :param polarities: Signs of recording poles
        """

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


    def compute_overall_CAP(self):
        """Simply sum over single fibre contributions. Result saved in ``self.CAP``.

        """

        arrayd_axonwise_CAPs = np.array(self.CAP_axonwise)

        # TODO: for variable time step, adjust here to minimum time step and then add.

        self.CAP = np.sum(arrayd_axonwise_CAPs, axis=0)


    def compute_single_axon_CAP(self, axon):
        """
        1. Calculate extracellular potential (LFP) for every electrode point defined in ``self.electrodePositions``
        2. Sum over points of single electrodes

        :param axon: The axon for which the single fibre action potential is computed.

        """

        # The contactpoints that constitute one electrode contact (e.g. a cuff electrode ring) have to be recovered,
        # summed up together per recording location along the axon

        # data structure: axon.imem[0] = t ; axon.imem[1] = imem Matrix

        # imemMatrix = axon.imem[1]
        # timeVector = axon.imem[0]

        # todo: changed this. Problem?
        CAP_axonwise = np.zeros(np.shape(axon.imem)[1])

        LFPmeans = []
        # colors = ['g', 'b']
        for poleIndex in range(self.numberOfPoles):

            # calculate LFP with LFPy from membrane currents for first pole
            # LFP = self.extPotMech.calculate_LFP(axon, self.electrodePositions[:, :, poleIndex])

            LFP = self.extPotMech.calculate_extracellular_potential(np.vstack([axon.xmid, axon.ymid, axon.zmid]).T, axon.imem, self.electrodePositions[:, :, poleIndex])
            CAP_axonwise += np.mean(LFP, 0)*self.polarities[poleIndex]

            # import matplotlib.pyplot as plt
            # plt.plot(np.mean(LFP,0))
            # plt.show()

        #     plt.plot(np.transpose(LFP), color=colors[poleIndex])
        #
        # plt.show()

        self.CAP_axonwise.append(CAP_axonwise)

    def clean_up(self):
        self.extPotMech = None



