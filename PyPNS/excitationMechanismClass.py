from abc import abstractmethod # ABCMeta,
import constants

class ExcitationMechanism(object):
    # __metaclass__ = ABCMeta

    def __init__(self):
        self.timeRes = constants.timeResStim

    @abstractmethod
    def connect_axon(self, axon):
        pass

    @abstractmethod
    def delete_neuron_objects(self):
        pass
