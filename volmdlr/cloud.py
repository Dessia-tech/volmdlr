# -*- coding: utf-8 -*-
"""

"""

import math
import numpy as npy
# import open3d

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

import volmdlr.primitives3d as p3d
from typing import List

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
        
    
    def to_shell(self, resolution: int = 10, normal = None, offset: float = 0):
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
        vec1, vec2 = xyz_vect[posmax-2], xyz_vect[posmax-1]
        subcloud2d_tosimp = [subcloud3d[n].to_2d(position_plane[n]*normal, vec1, vec2) for n in range(resolution)]
        subcloud2d = [sub.simplify(resolution=5) for sub in subcloud2d_tosimp]
        initial_polygon2d = [cloud2d.to_polygon() for cloud2d in subcloud2d]
        
        #Offsetting
        if offset != 0 :
            position_plane, initial_polygon2d = self.offset_to_shell(position_plane, initial_polygon2d, offset)
        
        areas = [0]*len(initial_polygon2d)
        for n, poly in enumerate(initial_polygon2d):
            if poly is not None :
                areas[n] = poly.area()
        avg_area = sum(areas)/len(areas)        
        
        polygon2d, polygon3d = [], []
        banned = []
        for n, poly in enumerate(initial_polygon2d):
            if poly is None or (poly.area()<avg_area/10 and (n not in [0,len(initial_polygon2d)-1])):
                resolution -= 1
                banned.append(n)
            else :
                polygon2d.append(poly)
                new_polygon = poly.to_3d(position_plane[n]*normal, vec1, vec2)
                polygon3d.append(new_polygon)
        [position_plane.pop(k) for k in banned[::-1]]
        faces = []
        
        for n in range(resolution):
            # print('sewing polygon', round(n/resolution*100, 2), '%')
            poly1 = polygon3d[n]
            poly1 = poly1.simplify()
            
            if n == resolution-1 or n == 0:
                plane3d = vmf.Plane3D.from_plane_vectors(position_plane[n]*normal, vec1, vec2)
                surf2d = vmf.Surface2D(polygon2d[n],[])
                surf2d = vmf.Surface2D(poly1.to_2d(position_plane[n]*normal, vec1, vec2),[])
                faces.append(vmf.PlaneFace3D(plane3d, surf2d))
            if n != resolution-1:
                poly2 = polygon3d[n+1]
                poly2 = poly2.simplify()
                
                coords = poly1.sewing(poly2, vec1, vec2)
                for trio in coords :
                    faces.append(vmf.Triangle3D(*trio))
        return vmf.ClosedShell3D(faces)
    
    # def alpha_shape(self, alpha:float, number_point_samples:int):
    #     '''
    #     Parameters
    #     ----------
    #     alpha : float
    #         the parameter alpha determines how precise the object surface reconstruction is wanted to be.
    #         The bigger the value of alpha is, more convex the final object will be. If it is smaller,
    #         the algorigthm is able to find the concave parts of the object, giving a more precise object
    #         surface approximation
    #     number_point_samples : int
    #         denotes the number of points to be used from the point cloud to reconstruct the surface. It uses poisson disk sampling algorithm
    #
    #     Returns
    #     -------
    #     Returns a ClosedShell3D object
    #
    #     '''
    #
    #     points = [[p.x, p.y, p.z] for p in self.points]
    #     array = npy.array(points)
    #     points = open3d.cpu.pybind.utility.Vector3dVector(array)
    #     pcd = open3d.geometry.PointCloud()
    #     pcd.points = points
    #     mesh = open3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd, alpha)
    #     mesh.compute_vertex_normals()
    #     if number_point_samples != None:
    #         pcd = mesh.sample_points_poisson_disk(number_point_samples)
    #         # tetra_mesh, pt_map = open3d.geometry.TetraMesh.create_from_point_cloud(pcd)
    #         mesh = open3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(
    #         pcd, alpha,
    #         # tetra_mesh, pt_map
    #         )
    #         mesh.compute_vertex_normals()
    #     # open3d.visualization.draw_geometries([mesh], mesh_show_back_face=True)
    #     vertices = [volmdlr.Point3D(float(x), float(y), float(z)) for x, y, z in list(npy.asarray(mesh.vertices))]
    #     triangles = [vmf.Triangle3D(vertices[p1], vertices[p2], vertices[p3], color = (1, 0.1, 0.1), alpha = 0.6) for p1, p2, p3 in list(npy.asarray(mesh.triangles))]
    #
    #     return vmf.ClosedShell3D(triangles)
        
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
    
    def extended_cloud(self, distance_extended : float): 
        #it works if distance_extended >= 0
        spheres, extended_points = [], []
        for pt in self.points :
            extended_zone = vm.p3d.Sphere(pt, distance_extended)
            sphere_primitive = extended_zone.shell_faces[0]
            
            spheres.append(vmf.ClosedShell3D([sphere_primitive]))
            
            extended_points.extend(sphere_primitive.triangulation().points)
            
        for sphere in spheres:
            clean_extended_zone = []
            for point in extended_points :
                if not sphere.point_belongs(point):
                    clean_extended_zone.append(point)
            extended_points = clean_extended_zone
            
        return extended_points
    
    def offset_to_shell(self, positions_plane: List[vmf.Plane3D], 
                        polygons2D: List[vmw.ClosedPolygon2D], offset: float):
        
        origin_f, origin_l = positions_plane[0], positions_plane[-1]
        
        new_position_plane = [origin_f-offset] + positions_plane[1:-1] + [origin_l+offset]
        
        new_poly = [poly.offset(offset) for poly in polygons2D]
        
        return new_position_plane, new_poly
    
class PointCloud2D(dc.DessiaObject):
    def __init__(self, points, name: str=''):
        self.points = points
        self.name = name
        
    def plot(self, ax=None, color='k'):
        if ax is None :
            fig, ax = plt.subplots()
        for pt in self.points :
            pt.plot(ax=ax, color=color)
        return ax
    
    def to_polygon(self):
        if not self.points:
            return None
        polygon = vmw.ClosedPolygon2D.points_convex_hull(self.points)
        # polygon = vmw.ClosedPolygon2D.concave_hull(self.points, -0.5, 0.0005)
        # polygon = vmw.ClosedPolygon2D.convex_hull_points(self.points)
        if polygon is None or math.isclose(polygon.area(), 0, abs_tol = 1e-6) :
            return None
        else : 
            return polygon
        
    def simplify(self, resolution = 5):
        if not self.points : 
            return PointCloud2D(self.points, name=self.name + '_none')
        
        x_list, y_list = [pt.x for pt in self.points], [pt.y for pt in self.points]
        xmin, xmax = min(x_list), max(x_list)
        ymin, ymax = min(y_list), max(y_list)
        
        x_slide = [xmin + n*(xmax-xmin)/(resolution-1) for n in range(resolution)]
        y_slide = [ymin + n*(ymax-ymin)/(resolution-1) for n in range(resolution)]
        
        points = []
        for x1, x2 in zip(x_slide, x_slide[1:]+[x_slide[0]]) :
            for y1, y2 in zip(y_slide, y_slide[1:]+[y_slide[0]]):
                box_points = []
                for pt in self.points :
                    if pt.x >= x1 and pt.x <= x2 :
                        if pt.y >= y1 and pt.y <= y2 :
                            box_points.append(pt)
                points.append(box_points)
        
        polys = [vmw.ClosedPolygon2D.points_convex_hull(pts) for pts in points]
        clean_points = []
        for poly in polys :
            if poly is not None :
                clean_points += poly.points
        return PointCloud2D(clean_points, name=self.name + '_clean')