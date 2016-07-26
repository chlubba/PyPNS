#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import neuron


class Cell:
    def __init__(self):
        neuron.h.dt = 0.025
        self.tstop = 20
        self.v_init = -70
        # sections
        self.soma = neuron.h.Section(name='soma', cell=self)
        self.soma.L = 30
        self.soma.diam = 30
        self.soma.nseg = 2
        self.dend = neuron.h.Section(name='dend', cell=self)
        self.dend.L = 1000
        self.dend.diam = 2
        self.dend.nseg = 20
        self.dend.connect(self.soma, 1, 0)
        self.allseclist = self.create_seclist()
        # mechanisms
        for sec in self.allseclist:
            sec.insert('pas')
        for sec in self.allseclist:
            sec.insert('extracellular')
        # insert field on each segment, initialize vmem, imem-recorders
        self.insert_vext()
        self.set_vmem_recorders()
        self.set_imem_recorders()

    def simulate(self):
        neuron.h.finitialize(self.v_init)
        neuron.h.fcurrent()
        while neuron.h.t < self.tstop:
            neuron.h.fadvance()
        # collect vmem and imem
        self.vmem = np.array(self.vmem)
        self.imem = np.array(self.imem)

    def create_seclist(self):
        allseclist = neuron.h.SectionList()
        allseclist.append(sec=self.soma)
        allseclist.append(sec=self.dend)
        return allseclist

    def set_vmem_recorders(self):
        self.vmem = neuron.h.List()
        for sec in self.allseclist:
            for seg in sec:
                memvrec = neuron.h.Vector()
                memvrec.record(seg._ref_v, neuron.h.dt)
                self.vmem.append(memvrec)

    def set_imem_recorders(self):
        self.imem = neuron.h.List()
        for sec in self.allseclist:
            for seg in sec:
                memirec = neuron.h.Vector()
                memirec.record(seg._ref_i_membrane, neuron.h.dt)
                self.imem.append(memirec)

    def insert_vext(self):
        # create some list of extracellular potentials on each segment and time vector
        t_ext = neuron.h.Vector(np.arange(20 / neuron.h.dt + 1) * neuron.h.dt)
        v_ext = []
        # play v_ext into e_extracellular reference
        for sec in self.allseclist:
            for seg in sec:
                # v_ext.append(neuron.h.Vector(np.random.rand(len(t_ext) - 0.5)))
                v_ext.append(neuron.h.Vector(np.random.rand(len(t_ext))))
                v_ext[-1].play(seg._ref_e_extracellular, t_ext)
        self.v_ext = v_ext
        self.t_ext = t_ext


if __name__ == '__main__':
    cell = Cell()
    cell.simulate()

    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    eim = ax1.matshow(cell.v_ext, cmap='spectral')
    cb1 = fig.colorbar(eim, ax=ax1)
    cb1.set_label('v_ext')
    ax1.axis(ax1.axis('tight'))
    iim = ax2.matshow(cell.imem, cmap='spectral')
    cb2 = fig.colorbar(iim, ax=ax2)
    cb2.set_label('imem')
    ax2.axis(ax2.axis('tight'))
    vim = ax3.matshow(cell.vmem, cmap='spectral')
    ax3.axis(ax3.axis('tight'))
    cb3 = fig.colorbar(vim, ax=ax3)
    cb3.set_label('vmem')
    ax3.set_xlabel('tstep')
    plt.show()