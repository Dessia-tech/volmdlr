# -*- coding: utf-8 -*-
"""

"""

import volmdlr.cloud
import volmdlr.core
import volmdlr as vm
import volmdlr.wires as vmw
import volmdlr.faces as vmf
import volmdlr.edges as vme
import matplotlib.pyplot as plt



# cloud3d = volmdlr.cloud.PointCloud3D.from_stl('')
cloud3d = volmdlr.cloud.PointCloud3D.from_stl('')
faces = cloud3d.subdescription_2d()
volum = volmdlr.core.VolumeModel(faces)
volum.babylonjs()

# points_poly1 = [vm.Point3D(-0.5268311157226563, -0.33678421020507815, 0.44895710584852433),
#                 vm.Point3D(-0.5354495239257813, -0.3411849060058594, 0.44895710584852433),
#                 vm.Point3D(-0.5429874877929688, -0.35194821166992185, 0.44895710584852433),
#                 vm.Point3D(-0.5456519165039062, -0.3583517150878906, 0.44895710584852433),
#                 vm.Point3D(-0.5558966064453125, -0.4118817138671875, 0.44895710584852433),
#                 vm.Point3D(-0.5547901000976563, -0.4208970947265625, 0.44895710584852433),
#                 vm.Point3D(-0.5523814086914063, -0.431831787109375, 0.44895710584852433),
#                 vm.Point3D(-0.54845458984375, -0.4380675048828125, 0.44895710584852433),
#                 vm.Point3D(-0.543478271484375, -0.44351849365234375, 0.44895710584852433),
#                 vm.Point3D(-0.5387916259765625, -0.44699050903320314, 0.44895710584852433),
#                 vm.Point3D(-0.5327277221679687, -0.4511383972167969, 0.44895710584852433),
#                 vm.Point3D(-0.528389404296875, -0.4534826049804688, 0.44895710584852433),
#                 vm.Point3D(-0.523892578125, -0.4556205139160156, 0.44895710584852433),
#                 vm.Point3D(-0.4628594970703125, -0.4817019958496094, 0.44895710584852433),
#                 vm.Point3D(-0.4197182006835937, -0.49059170532226565, 0.44895710584852433),
#                 vm.Point3D(-0.4072258911132812, -0.492669189453125, 0.44895710584852433),
#                 vm.Point3D(-0.39459478759765626, -0.4938609008789063, 0.44895710584852433),
#                 vm.Point3D(-0.38390960693359377, -0.494568603515625, 0.44895710584852433),
#                 vm.Point3D(-0.3722589111328125, -0.4946005859375, 0.44895710584852433),
#                 vm.Point3D(-0.360647705078125, -0.49353759765625, 0.44895710584852433),
#                 vm.Point3D(-0.3540591125488281, -0.49246218872070313, 0.44895710584852433),
#                 vm.Point3D(-0.34753570556640623, -0.4910741882324219, 0.44895710584852433),
#                 vm.Point3D(-0.22399929809570313, -0.46332778930664065, 0.4489571058485243),
#                 vm.Point3D(-0.22273770141601562, -0.46255389404296876, 0.44895710584852433),
#                 vm.Point3D(-0.21750790405273437, -0.45914111328125, 0.44895710584852433),
#                 vm.Point3D(-0.21426809692382812, -0.4562909851074219, 0.44895710584852433),
#                 vm.Point3D(-0.21527389526367188, -0.4526507873535156, 0.44895710584852433),
#                 vm.Point3D(-0.21610369873046875, -0.4500474853515625, 0.44895710584852433),
#                 vm.Point3D(-0.21691920471191406, -0.447577392578125, 0.44895710584852433),
#                 vm.Point3D(-0.21751730346679687, -0.4463346862792969, 0.44895710584852433),
#                 vm.Point3D(-0.218460205078125, -0.4452441101074219, 0.44895710584852433),
#                 vm.Point3D(-0.5189561157226562, -0.33894891357421875, 0.44895710584852433),
#                 vm.Point3D(-0.5228344116210938, -0.3377157897949219, 0.44895710584852433)]

# polygon1 = vmw.ClosedPolygon3D(points_poly1)

# points_poly2 = [vm.Point3D(-0.3888247985839844, -0.43719049072265626, 0.5092484079996744),
#                 vm.Point3D(-0.3965710144042969, -0.437202392578125, 0.5092484079996744),
#                 vm.Point3D(-0.41918109130859377, -0.43747021484375, 0.5092484079996744),
#                 vm.Point3D(-0.4192127990722656, -0.44862008666992187, 0.5092484079996744),
#                 vm.Point3D(-0.4155188903808594, -0.4788276062011719, 0.5092484079996744),
#                 vm.Point3D(-0.413130615234375, -0.47883810424804685, 0.5092484079996744),
#                 vm.Point3D(-0.4122532043457031, -0.47884100341796876, 0.5092484079996744),
#                 vm.Point3D(-0.4106119079589844, -0.4788460083007812, 0.5092484079996744),
#                 vm.Point3D(-0.408456787109375, -0.478852294921875, 0.5092484079996744),
#                 vm.Point3D(-0.3270098876953125, -0.47908309936523436, 0.5092484079996744),
#                 vm.Point3D(-0.18856179809570311, -0.47478271484375, 0.5092484079996744),
#                 vm.Point3D(-0.18802830505371093, -0.47472088623046876, 0.5092484079996744),
#                 vm.Point3D(-0.18574530029296876, -0.47419891357421873, 0.5092484079996744),
#                 vm.Point3D(-0.1839293975830078, -0.4729540100097656, 0.5092484079996744),
#                 vm.Point3D(-0.18271549987792968, -0.47109759521484373, 0.5092484079996744),
#                 vm.Point3D(-0.1822884979248047, -0.4688629150390625, 0.5092484079996744),
#                 vm.Point3D(-0.18222129821777344, -0.44926251220703123, 0.5092484079996744),
#                 vm.Point3D(-0.18262899780273437, -0.4470963134765625, 0.5092484079996744),
#                 vm.Point3D(-0.18382730102539063, -0.44523580932617185, 0.5092484079996744),
#                 vm.Point3D(-0.18563299560546875, -0.4439827880859375, 0.5092484079996744),
#                 vm.Point3D(-0.18774330139160156, -0.4435385131835938, 0.5092484079996744),
#                 vm.Point3D(-0.18793490600585938, -0.4435209045410156, 0.5092484079996744)]

# polygon2 = vmw.ClosedPolygon3D(points_poly2)

# vec1, vec2 = vm.X3D, vm.Y3D
# normal = vm.Z3D

# resolution = 20 #max([len(points_poly1),len(points_poly2)])
# primitives1 = [vme.LineSegment3D(point1, point2) for point1, point2 in 
#                zip(polygon1.points+[polygon1.points[0]], 
#                    polygon1.points[1:]+polygon1.points[:2])]
# primitives2 = [vme.LineSegment3D(point1, point2) for point1, point2 in 
#                zip(polygon2.points+[polygon2.points[0]], 
#                    polygon2.points[1:]+polygon2.points[:2])]
# contour1, contour2 = vmw.Contour3D(primitives1), vmw.Contour3D(primitives2)
# new_point1 = [contour1.point_at_abscissa(contour1.length()*n/(resolution-1)) for n in range(resolution)]
# new_point2 = [contour2.point_at_abscissa(contour2.length()*n/(resolution-1)) for n in range(resolution)]

# new_poly1, new_poly2 = vmw.ClosedPolygon3D(new_point1[:-1]), vmw.ClosedPolygon3D(new_point2[:-1])

# triangles = []
# for point1, point2, other_point in zip(new_poly1.points+[new_poly1.points[0]], 
#                                        new_poly1.points[1:]+new_poly1.points[:2],
#                                        new_poly2.points+[new_poly2.points[0]]):
#     triangles.append(vmf.Triangle3D(point1, point2, other_point))
        
# volum = volmdlr.core.VolumeModel(triangles)
# volum.babylonjs()
                         
# for point1, point2, other_point in zip(new_poly2.points+[new_poly2.points[0]],
#                                        new_poly2.points[1:]+new_poly2.points[:2],
#                                        new_poly1.points[1:]+new_poly1.points[:1]):
#     triangles.append(vmf.Triangle3D(point1, point2, other_point))
    
# volum = volmdlr.core.VolumeModel(triangles)
# volum.babylonjs()

# ax = contour1.plot()
# for pt in new_point1 :
#     pt.plot(ax=ax)
# new_point1[0].plot(ax=ax, color='r')
# new_point1[1].plot(ax=ax, color='b')
# new_point1[2].plot(ax=ax, color='m')
# new_point1[3].plot(ax=ax, color='g')
# contour2.plot(ax=ax)
# for pt in new_point2 :
#     pt.plot(ax=ax)
# new_point2[0].plot(ax=ax, color='r')
# new_point2[1].plot(ax=ax, color='b')
# new_point2[2].plot(ax=ax, color='m')
# new_point2[3].plot(ax=ax, color='g')    
