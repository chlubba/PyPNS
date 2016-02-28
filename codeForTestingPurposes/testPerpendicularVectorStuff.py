from math import pi, tan, cos, sin
import numpy as np
import scipy as sp

def perpendicular_vector(v):
    if v[0] == 0 and v[1] == 0:
        if v[2] == 0:
            raise ValueError('zero vector')
        # v is Vector(0, 0, v.z)
        v1 = np.array([0, 1, 0])

    v1 = np.array([-v[1], v[0], 0])
    v2 = np.cross(v,v1)

    randomAngle = np.random.uniform(0,2*pi,1)
    vRand = v1*cos(randomAngle) + v2*sin(randomAngle)

    return vRand/np.linalg.norm(vRand)


vector = np.array([1,0,0])
vP = perpendicular_vector(vector)
print vP
print np.linalg.norm(vP)
print np.inner(vP,vector)