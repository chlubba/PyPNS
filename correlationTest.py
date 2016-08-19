import numpy as np
import matplotlib.pyplot as plt

x = np.arange(0,10,0.1)
y1 = np.sin(x)
y2 = np.sin(x+np.pi/2)+10

corr = np.correlate(y1, y2, 'full')

plt.plot(corr[len(y1):])
plt.plot(y1)
plt.plot(y2)
plt.show()