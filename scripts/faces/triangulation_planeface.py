# -*- coding: utf-8 -*-
"""

"""

import volmdlr as vm
import volmdlr.faces as vmf
import volmdlr.stl as vmstl

list_tri = []
for stl_file in [
                'face0.stl',
                'face1.stl'
                  ]:
    file_faces = vmstl.Stl.from_file(stl_file)
    triangles = file_faces.triangles
    list_tri.append(triangles)
    # for face in file_faces.triangles:
    #     points = face.triangulation().points
    #     for three_pos in face.triangulation().triangles :
    #         triangles.append(vmf.Triangle3D(points[three_pos[0]], 
    #                                         points[three_pos[1]],
    #                                         points[three_pos[2]]))
  
for triangles in list_tri:
    print()          
    for tri in triangles : 
        print(tri.area())