import timeit

setup = '''
import numpy as np

numSteps = 10000.

x = [1.2191892,321.91892,50210378]
'''

print timeit.Timer('np.sqrt(np.sum(np.square(x)))', setup=setup).timeit(number=100)
print timeit.Timer('np.linalg.norm(x)', setup=setup).timeit(number=100)


