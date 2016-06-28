import numpy as np

field = np.array([[0., 0., -1.],
                  [0., 0., 0.],
                  [0., 0., 1.],
                  [0., 0., 2.],
                  [0., 1., -1.],
                  [0., 1., 0.],
                  [0., 1., 1.],
                  [0., 1., 2.],
                  [1., 0., -1.],
                  [1., 0., 0.],
                  [1., 0., 1.],
                  [1., 0., 2.],
                  [1., 1., -1.],
                  [1., 1., 0.],
                  [1., 1., 1.],
                  [1., 1., 2.],
                  ])

x = field[:,0]
y = field[:,1]
z = field[:,2]

# sort by coordinate values, x changing fastest, z slowest
orderIndices = np.lexsort((x, y, z))
x = x[orderIndices]
y = y[orderIndices]
z = z[orderIndices]

sortedField = np.vstack([x, y, z])

print sortedField