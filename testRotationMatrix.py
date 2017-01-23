import numpy as np

# def rotationMatrixFromVectors(a,b):
#     """compute rotation matrix that alignes b with a"""
#
#     lenA = np.linalg.norm(a)
#     if not lenA == 1:
#         a /= lenA
#
#     lenB = np.linalg.norm(b)
#     if not lenB == 1:
#         b /= lenB
#
#     v = np.cross(a,b)
#     s = np.linalg.norm(v)
#     c = np.dot(a,b)
#
#     vX = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
#
#     R = np.diag(np.ones(3)) + vX + np.dot(vX, vX)*1./(1. + c)

def rotationMatrixFromVectors(a,b):

    """np.dot(R,a) = b"""

    G = np.array([[np.dot(a, b), -np.linalg.norm(np.cross(a, b)), 0],
        [np.linalg.norm(np.cross(a, b)), np.dot(a, b), 0],
        [0, 0, 1]])
    F = np.array([a, (b - np.multiply(np.dot(a, b), a)) / np.linalg.norm(b - np.multiply(np.dot(a, b), a)), np.cross(b, a)])

    R = F.dot(G).dot(np.linalg.inv(F))

    return R

a = [1, 0, 0]
b = [0, 1, 0]

R = rotationMatrixFromVectors(a, b)

print a
print b
print np.linalg.norm(R)
print np.dot(R, b)