import numpy as np
import matplotlib.pyplot as plt


def func1(x, a, b, c, d):
    return a * (x + b) ** (-c) + d

z = np.linspace(0,0.003)
v = func1(z, 10**(-12), 0.0005, 2, 0)

plt.plot(z, v)
plt.show()