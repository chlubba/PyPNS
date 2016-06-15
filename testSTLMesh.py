# from stl import mesh
# from mpl_toolkits import mplot3d
# from matplotlib import pyplot
#
# # Create a new plot
# figure = pyplot.figure()
# axes = mplot3d.Axes3D(figure)
#
# # Load the STL files and add the vectors to the plot
# your_mesh = mesh.Mesh.from_file('/home/carl/Desktop/HalfDonut.stl')
# axes.add_collection3d(mplot3d.art3d.Poly3DCollection(your_mesh.vectors))
#
# # Auto scale to the mesh size
# scale = your_mesh.points.flatten(-1)
# axes.auto_scale_xyz(scale, scale, scale)
#
# # Show the plot to the screen
# pyplot.show()

import numpy as np
from stl import mesh

# Define the 8 vertices of the cube
vertices = np.array([\
    [-1, -1, -1],
    [+1, -1, -1],
    [+1, +1, -1],
    [-1, +1, -1],
    [-1, -1, +1],
    [+1, -1, +1],
    [+1, +1, +1],
    [-1, +1, +1]])
# Define the 12 triangles composing the cube
faces = np.array([\
    [0,3,1],
    [1,3,2],
    [0,4,7],
    [0,7,3],
    [4,5,6],
    [4,6,7],
    [5,1,2],
    [5,2,6],
    [2,3,6],
    [3,7,6],
    [0,1,5],
    [0,5,4]])

# vertices= [ [1, 0, 0],
#             [0, 1, 0],
#             [0, 0, 1]]

ringPoints = 10
angles = np.linspace(0, 2*np.pi, ringPoints, endpoint=False)

points = np.transpose(np.array(np.hstack([np.vstack([np.zeros(ringPoints),np.sin(angles), np.cos(angles)]),
                                         np.vstack([np.ones(ringPoints),np.sin(angles), np.cos(angles)])])))

# Create the mesh
cube = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        cube.vectors[i][j] = vertices[f[j],:]

# cube = mesh.Mesh(np.zeros(points.shape[0], dtype=mesh.Mesh.dtype))
# for i, f in enumerate(points):
#     cube.vectors[i] = points[i,:]

from mpl_toolkits import mplot3d
from matplotlib import pyplot

# Create a new plot
figure = pyplot.figure()
axes = mplot3d.Axes3D(figure)

# add the vectors to the plot
axes.add_collection3d(mplot3d.art3d.Poly3DCollection(cube.vectors))

# Auto scale to the mesh size
scale = cube.points.flatten(-1)
axes.auto_scale_xyz(scale, scale, scale)

# Show the plot to the screen
pyplot.show()