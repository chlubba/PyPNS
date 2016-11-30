import numpy as np

L = 0.01 # m
r1 = 200*10**-6 # m
r2 = 220*10**-6 # m
rho1 = 1/0.5
rho2 = 1/0.1

R = L / np.pi * (rho1 / (r1 ** 2) + rho2 / (r2 ** 2 - r1 ** 2))
print 'resistance = ' + str(R) + ' Ohm'


capacitances  = 10.**(-np.arange(6,13))

for c in capacitances:
    print 'c = ' + str(c) + ',tau = ' + str(R*c) + 's'

# ------------------ cap of wires

L = 0.5 # m
h = 1*10**-3 # distance in m
b = 0.5*10**-3 # radius of the wire

epsilon = 8.854*10**-12

C = np.pi*epsilon*L/(np.log(h/b))
print 'wire capacitance = ' + str(C*10**6) + ' uF'
