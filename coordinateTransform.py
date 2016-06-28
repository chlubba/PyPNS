import numpy as np

x = np.sin(np.arange(10)/10.*2*np.pi)
y = np.cos(np.arange(10)/10.*2*np.pi)
z = range(10)

rho = np.sqrt(x**2 + y**2)
phi = np.arctan2(y,x)

print rho
print phi

# --------------------- timer ----------------------------------------------

import timeit

setup = '''
import numpy as np

numSteps = 10000.

x = np.sin(np.arange(numSteps)/numSteps*2*np.pi)
y = np.cos(np.arange(numSteps)/numSteps*2*np.pi)
'''

print timeit.Timer('np.sqrt(x**2 + y**2)', setup=setup).timeit(number=100)
print timeit.Timer('np.arctan2(y,x)', setup=setup).timeit()

