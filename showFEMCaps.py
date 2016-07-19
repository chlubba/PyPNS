import numpy as np
import matplotlib.pyplot as plt

filenames = ['/media/carl/4ECC-1C44/PyPN/FEM_CAPs/z0.05m_100steps_2.npy',
             '/media/carl/4ECC-1C44/PyPN/FEM_CAPs/z0.003m_100steps_2.npy',
             ] # '/media/carl/4ECC-1C44/PyPN/FEM_CAPs/z0.001m_500steps_2.npy',
             # '/media/carl/4ECC-1C44/PyPN/FEM_CAPs/z0.001m_50000steps.npy',

labels = ['500 $\mu m$', '30 $\mu m$'] # , '2 $\mu m$'] # , '0.02 $\mu m$'

for index, filename in enumerate(filenames):
    electrodeData = np.load(filename)
    plt.plot(np.sum(electrodeData, axis=0), label=labels[index])

plt.legend()
plt.show()