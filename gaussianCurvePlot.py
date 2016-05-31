import numpy as np
import matplotlib.pyplot as plt

sigma = 0.3
mu = 0.7

sigma = 0.6
mu = 2.0

s = np.random.normal(mu, sigma, 1000)

count, bins, ignored = plt.hist(s, 10, normed=True)
plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu)**2 / (2 * sigma**2) ),linewidth=2, color='r')
# plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6])
plt.xticks(np.arange(0,4,0.2))
plt.show()