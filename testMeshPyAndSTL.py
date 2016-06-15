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

