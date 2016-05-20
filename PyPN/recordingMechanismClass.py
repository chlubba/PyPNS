# from abc import ABCMeta, abstractmethod

class RecordingMechanism(object):
    # __metaclass__ = ABCMeta

    def __init__(self, numberOfElectrodes):

        self.numberOfElectrodes = numberOfElectrodes

    # @abstractmethod
    def setup_recording_elec(self):
        pass

    # @abstractmethod
    def set_number_of_elecs(self):
        pass

    # def (self):
    #     pass

class CuffElectrodes(RecordingMechanism):

    def __init__(self, ):

        numberOfElectrodes = []

        super(CuffElectrodes,self).__init__(numberOfElectrodes)

    def setup_recording_elec(self):
        pass



elec = CuffElectrodes()

print elec