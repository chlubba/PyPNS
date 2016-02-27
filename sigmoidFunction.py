import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


x = np.arange(0,1,0.01)
y = 1-1/(1+np.exp(-20*(x - 0.5)))

plt.plot(x,y)
plt.show()