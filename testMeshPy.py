import meshpy
import meshpy.geometry
import numpy as np

points, facets, holelists, markers = meshpy.geometry.make_cylinder(radius=10, height=10, radial_subdivisions=10, height_subdivisions=1)
# # print cylinder
from meshpy.tet import MeshInfo, build
# import sfepy.discrete.fem.mesh

mesh_info = MeshInfo()
mesh_info.set_points(points)
mesh_info.set_facets(np.squeeze(np.array(facets)))

mesh = build(mesh_info)

# mesh_info = MeshInfo()
# mesh_info.set_points([
#     (0,0,0), (2,0,0), (2,2,0), (0,2,0),
#     (0,0,12), (2,0,12), (2,2,12), (0,2,12),
#     ])
# mesh_info.set_facets([
#     [0,1,2,3],
#     [4,5,6,7],
#     [0,4,5,1],
#     [1,5,6,2],
#     [2,6,7,3],
#     [3,7,4,0],
#     ])
# mesh = build(mesh_info)

print "Mesh Points:"
for i, p in enumerate(mesh.points):
    print i, p
print "Point numbers in tetrahedra:"
for i, t in enumerate(mesh.elements):
    print i, t
mesh.write_vtk("test.vtk")

from stl import mesh as STLMesh

vectors = np.zeros((mesh.elements.__len__(), 4, 3))

# Create the mesh
cylinder = STLMesh.Mesh(np.zeros(mesh.elements.__len__(), dtype=STLMesh.Mesh.dtype))
for i, f in enumerate(mesh.elements):
    for j in range(4):
        # cylinder.vectors[i][j] = mesh.points[f[j]]
        vectors[i][j] = mesh.points[f[j]]

cylinder.vectors = vectors

# Write the mesh to file "cube.stl"
cylinder.save('cylinder.stl')


from mpl_toolkits import mplot3d
from matplotlib import pyplot

# Create a new plot
figure = pyplot.figure()
axes = mplot3d.Axes3D(figure)

# add the vectors to the plot
axes.add_collection3d(mplot3d.art3d.Poly3DCollection(cylinder.vectors))

# Auto scale to the mesh size
scale = cylinder.points.flatten(-1)
axes.auto_scale_xyz(scale, scale, scale)

# Show the plot to the screen
pyplot.show()

