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


p1 = vm.Point2D([1,1])
p2 = vm.Point2D([2,1])
p3 = vm.Point2D([2.5,2])
p4 = vm.Point2D([1.5,3])
p5 =  vm.Point2D([1,2])
p6 =  vm.Point2D([1.6666666666666667, 2.3333333333333335])


#tr1=vmmesh.TriangularElement([p3,p4,p6])
#tr1.plot()

l1 = vm.LineSegment2D(p1,p2)
l2 = vm.LineSegment2D(p2,p3)
l3=vm.LineSegment2D(p3,p4)
l4 = vm.LineSegment2D(p4,p5) 
l5 = vm.LineSegment2D(p5,p1)
contour=vm.Contour2D([l1,l2,l3,l4,l5])
polygon=contour.polygon
mesher=vmmesh.Mesher(contour,[],[],1.5 )
triangles=mesher.triangulation_polygone_recursive(polygon)
all_triangles=mesher.assemble_mesh(triangles,0.001)




print(len(all_triangles))
print(len(set(all_triangles)))
