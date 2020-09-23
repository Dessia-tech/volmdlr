import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as npy
import matplotlib.tri as mtri
import volmdlr as vm
import volmdlr.mesh as vmmesh
import finite_elements.elasticity as els
import finite_elements.core as corefe
import math
from scipy import sparse
from scipy import linalg
from finite_elements.core import steel,aluminium
import time 
from dessia_common import DessiaObject
from typing import TypeVar, List, Tuple
from scipy.spatial import Delaunay
from itertools import product 
p1 = vm.Point2D([1,1])
p2 = vm.Point2D([2,1])
p3 = vm.Point2D([2.5,2])
p4 = vm.Point2D([1.5,3])
p5 =  vm.Point2D([1,2])
p6 =  vm.Point2D([1.6666666666666667, 2.3333333333333335])
points=npy.array([[1,1],[2,1],[2.5,2],[1.5,3],[1,2]])

tr1=vmmesh.TriangularElement([p3,p4,p6])



l1 = vm.LineSegment2D(p1,p2)
l2 = vm.LineSegment2D(p2,p3)
l3=vm.LineSegment2D(p3,p4)
l4 = vm.LineSegment2D(p4,p5) 
l5 = vm.LineSegment2D(p5,p1)
contour=vm.Contour2D([l1,l2,l3,l4,l5])

# L=[p1,p2]
# M=[p3,p4]
# print(list(product(L,M))[0][0])



polygon=contour.polygon
all_points=polygon.points

mesher=vmmesh.Mesher(contour,[],1.7)
polygons=polygon.delaunay_triangulation()
triangles=mesher.polygon_to_triangles(polygons)



all_triangles=mesher.assemble_mesh(triangles)

print(len(all_triangles))
print(len(set(all_triangles)))


# all_triangles=mesher.assemble_mesh(triangles)

element_group=vmmesh.ElementsGroup(all_triangles,'element_group')

mesh=vmmesh.Mesh([element_group])
mesh.plot()  





                                              
#all_triangles=mesher.assemble_mesh(triangles,0.1)
# element_group=vmmesh.ElementsGroup(all_triangles,'element_group')
# mesh=vmmesh.Mesh([element_group])
# mesh.plot()  
# print(len(all_triangles))
# print(len(set(all_triangles)))





