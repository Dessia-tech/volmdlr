# -*- coding: utf-8 -*-
"""

"""

import math
import numpy as npy


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import volmdlr as vm
import volmdlr.core
import volmdlr.wires as vmw
import volmdlr.faces as vmf
import volmdlr.step as vstep
import volmdlr.stl as vmstl
import volmdlr.edges as vme
import dessia_common as dc

class PointCloud3D(dc.DessiaObject):
    def __init__(self, points, name: str=''):
        self.points = points
        self.name = name
        
    @classmethod
    def from_stl(cls, file_path):
        list_points = vmstl.Stl.from_file(file_path).extract_points_BIS()
        
        return cls(list_points, name='from_stl')
    
    def _bounding_box(self):
        return vm.core.BoundingBox.from_points(self.points)
    
    def to_2d(self, plane_origin, x, y):
        list_points2d = [pt3d.to_2d(plane_origin, x, y) for pt3d in self.points]
        return PointCloud2D(list_points2d, name='3d_to_2d')
    
    def extract(self, u, umin, umax):# -> List[PointCloud3D] :
        extracted_points = []
        for points in self.points :
            dist_to_plane = points.dot(u)
            if dist_to_plane > umin and dist_to_plane < umax :
                extracted_points.append(points)
        return PointCloud3D(extracted_points)
        
    
    def to_shell(self, resolution = 10, normal = None):
        #normal has to be a fondamental vector : X3D, Y3D or Z3D
        bbox = self._bounding_box()
        xyz_bbox = [[bbox.xmin, bbox.xmax], [bbox.ymin,bbox.ymax], [bbox.zmin,bbox.zmax]]
        xyz_list = [l[1]-l[0] for l in xyz_bbox]
        absxyz_list, xyz_vect = [abs(length) for length in xyz_list], [vm.X3D, vm.Y3D, vm.Z3D]

        if normal is None :
            posmax = xyz_list.index(max(absxyz_list))
            normal = xyz_vect[posmax]
            
        else :
            posmax = 0
            for n, vect in enumerate(xyz_vect):
                if vect == normal :
                    posmax = n
        dist_between_plane = xyz_list[posmax]/(resolution-1)
        position_plane = [xyz_bbox[posmax][0] + n*dist_between_plane for n in range(resolution)]
        
        subcloud3d = [self.extract(normal, pos_plane-dist_between_plane/2, pos_plane+dist_between_plane/2) for pos_plane in position_plane]
        print('subcloud3D CREATED')
        vec1, vec2 = xyz_vect[posmax-2], xyz_vect[posmax-1]
        subcloud2d = [subcloud3d[n].to_2d(position_plane[n]*normal, vec1, vec2) for n in range(resolution)]
        # return subcloud2d
        # fig = plt.figure()
        # ax = fig.add_subplot(111)
        # for point in subcloud2d[0].points[1::10]:
            # point.plot(ax=ax)
        print('subcloud2D CREATED')
        print('CREATING POLYGONS')
        initial_polygon2d = [cloud2d.to_polygon() for cloud2d in subcloud2d]
        
        polygon2d, polygon3d = [], []
        for pos_plane, poly in zip(position_plane, initial_polygon2d) :
            if poly is None :
                resolution -= 1
            else :
                polygon2d.append(poly)
                polygon3d.append(poly.to_3d(pos_plane*normal, vec1, vec2))
        # return polygon3d
                
        faces = []
        max_poly_resolution = int(sum([len(poly.points) for poly in polygon3d])/len(polygon3d))+1
        
        fig = plt.figure()
        ax = Axes3D(fig)
        for poly in polygon3d :
            poly.plot(ax=ax)
        
        for n in range(resolution):
            print('sewing polygon', round(n/resolution*100, 2), '%')
            poly1 = polygon3d[n]
            poly1 = poly1.simplify(0.05, 0.1)
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            for point in poly1.points:
                point.plot(ax= ax, color = 'g')
            for line in poly1.line_segments:
                line.plot(ax= ax)
            if n == resolution-1 or n == 0:
                plane3d = vmf.Plane3D.from_plane_vectors(position_plane[n]*normal, vec1, vec2)
                surf2d = vmf.Surface2D(polygon2d[n],[])
                faces.append(vmf.PlaneFace3D(plane3d, surf2d))
            if n != resolution-1:
                poly2 = polygon3d[n+1]
                poly2 = poly2.simplify(0.05, 0.1)
                # coords = poly1.sewin1g_with(poly2, vec1, vec2, normal, resolution = max_poly_resolution)
                coords = poly1.sewing(poly2)
                for trio in coords :
                    faces.append(vmf.Triangle3D(trio[0], trio[1], trio[2]))
                # volum = volmdlr.core.VolumeModel(faces)
                # volum.babylonjs()
                # fig = plt.figure()
                # ax = fig.add_subplot(111, projection='3d')
                # poly1.plot(ax=ax, color='g')
                # poly2.plot(ax= ax, color = 'r')
                # for point in poly1.points + poly2.points:
                #     point.plot(ax= ax, color = 'b')
                
        
        return vmf.ClosedShell3D(faces)

    @classmethod        
    def from_step(cls, step_file:str):
        step = vstep.Step(step_file)
        points = step.to_points()
        return cls(points)
    
    def plot(self, ax = None):
        ax = self.points[0].plot(ax = ax)
        for point in self.points[1::1000]:
            point.plot(ax = ax)
            
        return ax
    
class PointCloud2D(dc.DessiaObject):
    def __init__(self, points, name: str=''):
        self.points = points
        self.name = name
        
    def plot(self):
        fig, ax = plt.subplots()
        for pt in self.points :
            pt.plot(ax=ax)
        return ax
    
    def to_polygon(self):
        polygon = vmw.ClosedPolygon2D.points_convex_hull(self.points)
        # polygon = vmw.ClosedPolygon2D.concave_hull(self.points, 0.2, 0.000005)
        # polygon = vmw.ClosedPolygon2D.convex_hull_points(self.points)
        # fig = plt.figure()
        # ax = fig.add_subplot(111)
        # for point in polygon.points:
        #     point.plot(ax= ax, color = 'g')
        # for line in polygon.line_segments:
        #     line.plot(ax= ax)
        # polygon.plot()
        # print('polygon:', polygon)
        if polygon is None or math.isclose(polygon.area(), 0, abs_tol = 1e-6) :
            return None
        else : 
            return polygon