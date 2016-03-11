import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0,1,1000)

y = 1/(1+np.exp(-20*(x - 0.7)))

plt.plot(x,y)
plt.show()

