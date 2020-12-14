#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 10:29:52 2020

@author: mehdigasmi
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as npy
import matplotlib.tri as mtri
import volmdlr.edges as edges
import volmdlr.wires as wires
import volmdlr.faces as faces
import volmdlr.core_compiled as vm
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


p1=vm.Point2D(0,0)
p2=vm.Point2D(2,2)
p3=vm.Point2D(3,4)
p4=vm.Point2D(5,4)
p5=vm.Point2D(2.5,3.5)
p6=vm.Point2D(1,1.5)
p7=vm.Point2D(4,4.2)
p8=vm.Point2D(2.1,1.5)
p9=vm.Point2D(4,5)
p10=vm.Point2D(4,3.5)
l8=edges.LineSegment2D(p4, p1)
l1=edges.LineSegment2D(p1, p2)
a=edges.Arc2D(p2, p5, p3)
l2=edges.LineSegment2D(p3, p4)
w=wires.Wire2D([l1,a,l2])

a_2=edges.Arc2D(p1, p6, p2)
l3=edges.LineSegment2D(p2, p3)
l4=edges.LineSegment2D(p3, p4)
w2=wires.Wire2D([a_2,l3,l4])
w2.plot()
ax=w.plot()
offset=w.offset(-1)
offset.plot(ax=ax)
ax=w2.plot()
offset_2=w2.offset(-1)
offset_2.plot(ax=ax)
a_3=edges.Arc2D(p3, p7, p4)
l5=edges.LineSegment2D(p3, p4)
w3=wires.Wire2D([l1,l3,a_3])
w3.plot()
ax=w3.plot()
offset_3=w3.offset(-1)
offset_3.plot(ax=ax)

w4=wires.Wire2D([a_2,a,l2])
ax=w4.plot()
offset_4=w4.offset(-1)
offset_4.plot(ax=ax)

w5=wires.Wire2D([l1,a,a_3])
ax=w5.plot()
offset_5=w5.offset(-1)
offset_5.plot(ax=ax)

a_5=edges.Arc2D(p1,p8,p2)
w6=wires.Wire2D([a_5,l3,l2])
ax=w6.plot()
offset_6=w6.offset(-1)
offset_6.plot(ax=ax)

w7=wires.Contour2D([a_2,a,a_3,l8])
w_0=wires.Contour2D([l8,a,a_3,a_2])
print(w_0==w7)
ax=w7.plot()

offset_7=w7.offset(-0.2)
print(offset_7.primitives)
print(offset_7.self_intersects())
rep=offset_7.repair_single_intersection()
for p in rep:
    p.plot()
rep2=offset_7.repair_intersections([offset_7],[])
for p in rep2:
    p.plot()
print(rep2[2]==rep2[3])
rep2[2].plot()
rep2[3].plot()
offset_7.plot()






l6=edges.LineSegment2D(p3, p9)
w8=wires.Wire2D([a_2,a,l6])
ax=w8.plot()
offset_8=w8.offset(-1)
offset_8.plot(ax=ax)

l7=edges.LineSegment2D(p2, p1)
a_5=edges.Arc2D(p2, p5, p3)
l8=edges.LineSegment2D(p4, p3)
w9=wires.Wire2D([l7,a_5,l8])
ax=w9.plot()
offset_9=w9.offset(-1)
offset_9.plot(ax=ax)
ax=w9.plot()
for p in w9.ordered_primitives():
    p.start.plot(ax=ax,color='r')
    p.end.plot(ax=ax,color='b')
    
a_6=edges.Arc2D(p3, p10, p4)
w10=wires.Wire2D([l1,l3,a_6])
ax=w10.plot() 
offset_10=w10.offset(0.25)
offset_10.plot(ax=ax)
    
l1=edges.LineSegment2D(p1, p2)
a=edges.Arc2D(p2, p5, p3)
l2=edges.LineSegment2D(p3, p4)

w11=wires.Contour2D([l1,a,l2,l8]) 
ax=w11.plot() 
offset_11=w11.offset(-1.1)
offset_11.plot(ax=ax) 
    
# if len(intersections)==1:
#                offset_intersections.append(([intersections[0],intersections[0]],'Line',i))
#            else :
#                 if intersections[0].point_distance(self.primitives[i].interior)>intersections[1].point_distance(self.primitives[i].interior):
                                                                                                               
#                     intersections.reverse()
#                 offset_intersections.append((intersections,'Circles',i))  
                
#  if intersections[0].point_distance(self.primitives[i].end)>intersections[1].point_distance(self.primitives[i].end):
#                          intersections.reverse()
#                     offset_intersections.append((intersections,'Circle',i+1))
                    
# if intersections[0].point_distance(self.primitives[i+1].start)>intersections[1].point_distance(self.primitives[i+1].start):
#                         intersections.reverse()
#                      offset_intersections.append((intersections,'Line',i))                    