from abc import ABCMeta, abstractmethod

class ExcitationMechanism(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def connect_axon(self, axon):
        pass

    @abstractmethod
    def delete_neuron_objects(self):
        pass