import numpy as np

f = open("/home/carl/testSavetxt/test.dat", 'ab')

matrix = np.multiply(np.ones([5,5]),5)

np.savetxt(f, matrix)

f.close()