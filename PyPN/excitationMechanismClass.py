from abc import abstractmethod # ABCMeta,

class ExcitationMechanism(object):
    # __metaclass__ = ABCMeta

    def __init__(self):
        self.timeRes = 0

    @abstractmethod
    def connect_axon(self, axon):
        pass

    @abstractmethod
    def delete_neuron_objects(self):
        pass