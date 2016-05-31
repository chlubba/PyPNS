import numpy as np

a = np.array([]).reshape(0,8)

a = np.vstack([a, np.ones(8)])
a = np.vstack([a, np.ones(8)*2])

print a