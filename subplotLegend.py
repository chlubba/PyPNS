import matplotlib.pyplot as plt
import numpy as np

f, (ax1, ax2) = plt.subplots(2,1)

ax1.plot(np.arange(5), label='lala')
ax1.legend()

ax2.plot(range(5), label='lala')
ax2.plot(np.ones(5), label='lala')
ax2.legend()

plt.show()