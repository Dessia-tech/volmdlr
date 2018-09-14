#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 14:08:23 2017

@author: steven
"""
import numpy as npy

import volmdlr
import math


class RoundedLineSegments3D(volmdlr.Wire3D):
    def __init__(self, points, radius, closed=False, name=''):        
        self.points=points
        self.radius=radius
        self.closed=closed
        # Construncting Arcs and lines of profile
        if closed:
            p1s=[points[-1]]+points[:-1]
            pis=points
            p2s=points[1:]+[points[0]]
        else:
            p1s=[points[-1]]+points[:-2]
            pis=points[:-1]
            p2s=points[1:]
    
        points_l=points[:]
#        primitives=[]
        arcs=[]
        for i in range(len(points)):
            if i in radius:
                r=self.radius[i]
                p1=p1s[i]
                pi=pis[i]
                p2=p2s[i]

                dist1 = (p1 - pi).Norm()
                dist2 = (p2 - pi).Norm()
                dist3 = (p1 - p2).Norm()
                alpha = math.acos(-(dist3**2-dist1**2-dist2**2)/(2*dist1*dist2))/2.
                dist = r/math.tan(alpha)

                u1 = (p1 - pi) / dist1
                u2 = (p2 - pi) / dist2
                
                p3 = pi + u1*dist
                p4 = pi + u2*dist

                n = u1.Cross(u2)
                n /= n.Norm()
                v1 = u1.Cross(n)
                v2 = u2.Cross(n)
                
                l1 = volmdlr.Line3D(p3, p3+v1)
                l2 = volmdlr.Line3D(p4, p4+v2)
                c, _ = l1.MinimumDistancePoints(l2)
                
                u3 = u1 + u2# mean of v1 and v2
                u3 /= u3.Norm()

                interior = c - u3 * r 
    
                points_l[i] = (p3, p4)                           
                arcs.append(volmdlr.Arc3D(p3, interior, p4))


        lines=[]

        last_point=points_l[0]
        for p in points_l[1:]+[points_l[0]]:
            if type(p)==tuple:
                lines.append(volmdlr.LineSegment3D(last_point,p[0]))
                last_point=p[1]
            else:
                lines.append(volmdlr.LineSegment3D(last_point,p))
                last_point=p
        if not closed:
            del lines[-1]
            
        primitives=lines[:]
        nlines=len(lines)
        narcs=len(arcs)
        for ii,i in enumerate(sorted(radius.keys(),reverse=True)):
            if i>nlines:
                primitives.append(arcs[narcs-ii-1])
            else:
                primitives.insert(i,arcs[narcs-ii-1])
        
        volmdlr.CompositePrimitive3D.__init__(self, primitives, name)        
        
        
    def Rotation(self, center, angle, copy=True):
        if copy:
            return RoundedLineSegments3D([p.Rotation(center,angle,copy=True) for p in self.points],self.radius,self.closed,self.name)
        else:
            self.__init__([p.Rotation(center,angle,copy=True) for p in self.points],self.radius,self.closed,self.name)
            
    def Translation(self,offset,copy=True):
        if copy:
            return RoundedLineSegments3D([p.Translation(offset,copy=True) for p in self.points],self.radius,self.closed,self.name)
        else:
            self.__init__([p.Translation(offset,copy=True) for p in self.points],self.radius,self.closed,self.name)

    def MPLPlot(self, ax):
        for primitive in self.basis_primitives:
            primitive.MPLPlot(ax)
        
class Cylinder(volmdlr.Primitive3D):
    def __init__(self, position, axis, radius, width, name=''):
        volmdlr.Primitive3D.__init__(self, name)
        self.position = position
        axis.Normalize()
        self.axis = axis
        self.radius = radius
        self.width = width
        
    def Volume(self):
        return self.width * math.pi * self.radius**2
        
    def FreeCADExport(self, ip):
        if self.radius>0:
            name='primitive'+str(ip)
            e=str(1000*self.width)
            r=str(1000*self.radius)
            position=1000*(self.position - self.axis*self.width/2.)
            x,y,z=position
            x=str(x)
            y=str(y)
            z=str(z)
    
            ax,ay,az=self.axis
            ax=str(ax)
            ay=str(ay)
            az=str(az)
            return name+'=Part.makeCylinder('+r+','+e+',fc.Vector('+x+','+y+','+z+'),fc.Vector('+ax+','+ay+','+az+'),360)\n'
        else:
            return ''
    
    def Babylon(self):
        ya,xa,za=self.axis
        theta=math.acos(za/self.width)
        phi=math.atan(ya/xa)
        x,z,y=self.position
        s='var cylinder = BABYLON.Mesh.CreateCylinder("{}", {}, {}, {}, 30, 1, scene,false, BABYLON.Mesh.DEFAULTSIDE);'.format(self.name,self.width,2*self.radius,2*self.radius)
        s+='cylinder.position = new BABYLON.Vector3({},{},{});\n;'.format(x,y,z)
        s+='cylinder.rotation.x={}\n;'.format(-theta*math.sin(phi))
        s+='cylinder.rotation.y={}\n;'.format(theta*math.cos(phi))
        s+='cylinder.rotation.z={}\n;'.format(phi)
        return s
    
class HollowCylinder(volmdlr.Primitive3D):
    def __init__(self, position, axis, inner_radius, outer_radius, width, name=''):
        volmdlr.Primitive3D.__init__(self, name)
        self.position = position
        axis.Normalize()
        self.axis = axis
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.width = width
        
    def Volume(self):
        return self.width * math.pi* (self.outer_radius**2 - self.inner_radius**2)


    def FreeCADExport(self, ip):
        if self.outer_radius > 0.:
            name = 'primitive'+str(ip)
            re=str(1000*self.outer_radius)
            ri=str(1000*self.inner_radius)        
            position = self.position-self.axis*self.width/2
            x,y,z=1000*position
            ax,ay,az=self.axis
            x=str(x)
            y=str(y)
            z=str(z)
            ax=str(ax)
            ay=str(ay)
            az=str(az)
    
            s='C2= Part.makeCircle('+re+',fc.Vector('+x+','+y+','+z+'),fc.Vector('+ax+','+ay+','+az+'))\n'
            s+='W2=Part.Wire(C2.Edges)\n'
            s+='F2=Part.Face(W2)\n'
            
            if self.inner_radius!=0.:
                s+='C1= Part.makeCircle('+ri+',fc.Vector('+x+','+y+','+z+'),fc.Vector('+ax+','+ay+','+az+'))\n'        
                s+='W1=Part.Wire(C1.Edges)\n'
                s+='F1=Part.Face(W1)\n'        
                s+='F2=F2.cut(F1)\n'        
    
            vx,vy,vz=self.axis*self.width*1000
            vx=str(vx)
            vy=str(vy)
            vz=str(vz)
            
            s+=name+'=F2.extrude(fc.Vector('+vx+','+vy+','+vz+'))\n'
            return s
        
        else:
            return ''
    
    def Babylon(self):
        ya,xa,za=self.axis# to counter y definition in babylon
        theta=math.acos(za/self.width)
        phi=math.atan(ya/xa)
        x,z,y=self.position
        s='var cylinder = BABYLON.Mesh.CreateCylinder("{}", {}, {}, {}, 30, 1, scene,false, BABYLON.Mesh.DEFAULTSIDE);'.format(self.name,self.width,2*self.outer_radius,2*self.outer_radius)
        s+='cylinder.position = new BABYLON.Vector3({},{},{});\n;'.format(x,y,z)
        s+='cylinder.rotation.x={}\n;'.format(-theta*math.sin(phi))
        s+='cylinder.rotation.y={}\n;'.format(theta*math.cos(phi))
        s+='cylinder.rotation.z={}\n;'.format(phi)
        return s
    
class ExtrudedProfile(volmdlr.Primitive3D):
    """

    """
    def __init__(self, plane_origin, x, y, outer_contour2d, inner_contours2d,
                 extrusion_vector, name=''):
        volmdlr.Primitive3D.__init__(self, name)
        self.outer_contour2d = outer_contour2d
        self.outer_contour3d = outer_contour2d.To3D(plane_origin, x, y)
        self.inner_contours2d = inner_contours2d
        self.extrusion_vector = extrusion_vector
        self.inner_contours3d = []
        self.x = x
        self.y = y
        for contour in inner_contours2d:
            self.inner_contours3d.append(contour.To3D(plane_origin, x, y))
        
    def MPLPlot(self, ax):
        for contour in self.contours3D:
            for primitive in contour:
                primitive.MPLPlot(ax)
        
    def FreeCADExport(self, ip):
        name='primitive'+str(ip)
        
        s = 'Wo = []\n'
        s += 'Eo = []\n'
        for ip, primitive in enumerate(self.outer_contour3d.basis_primitives):
            s += primitive.FreeCADExport('L{}'.format(ip))
            s += 'Eo.append(Part.Edge(L{}))\n'.format(ip)
        s += 'Wo.append(Part.Wire(Eo[:]))\n'
        s += 'Fo = Part.Face(Wo)\n'
        
        s += 'Fi = []\n'
        s += 'W = []\n'
        for ic,contour in enumerate(self.inner_contours3d): 
            s+='E = []\n'
            for ip, primitive in enumerate(contour.basis_primitives):
                s += primitive.FreeCADExport('L{}_{}'.format(ic, ip))
                s += 'E.append(Part.Edge(L{}_{}))\n'.format(ic, ip)
            s += 'Wi = Part.Wire(E[:])\n'
            s += 'Fi.append(Part.Face(Wi))\n'
            
        if len(self.inner_contours3d) != 0:
            s += 'Fo = Fo.cut(Fi)\n'
        e1, e2, e3 = self.extrusion_vector
        e1=e1*1000
        e2=e2*1000
        e3=e3*1000
        
        s+='{} = Fo.extrude(fc.Vector({}, {}, {}))\n'.format(name,e1,e2,e3)
        return s
    
    def Area(self):
        areas=[c.Area() for c in self.contours2D]
        sic=list(npy.argsort(areas))[::-1]# sorted indices of contours
        area=areas[sic[0]]
        
        for i in sic[1:]:
            area-=self.contours2D[i].Area()
        return area
    
    def Volume(self):
#        e=extrusion_vector/norm(extrusion_vector)
        z = npy.cross(self.x.vector,self.y.vector)
        z.Normalize()
        coeff = npy.dot(self.extrusion_vector, z)
        
        return self.Area()*coeff
        
        

class Sphere(volmdlr.Primitive3D):
    def __init__(self,center,radius,name=''):
        volmdlr.Primitive3D.__init__(self,name)
        self.center = center
        self.radius = radius
        self.position = center
    
    def Volume(self):
        return 4/3*math.pi*self.radius**3
    
    def FreeCADExport(self, ip, ndigits=3):
        name = 'primitive'+str(ip)
        r = 1000*self.radius
        x, y, z = npy.round(1000*self.center.vector,ndigits)
        return '{}=Part.makeSphere({},fc.Vector({},{},{}))\n'.format(name,r,x,y,z)


class RevolvedProfile(volmdlr.Primitive3D):
    """
    
    """
    def __init__(self, plane_origin, x, y, contours2d, axis_point, 
                 axis,angle=2*math.pi, name=''):
        volmdlr.Primitive3D.__init__(self, name)
        self.contours2d = contours2d
        self.axis_point = axis_point
        self.axis = axis
        self.angle = angle
        self.plane_origin = plane_origin
        self.x = x
        self.y = y
        
        self.contours3d = []
        for contour in contours2d:
            self.contours3d.append(contour.To3D(plane_origin, x, y))
        
    def MPLPlot(self, ax):
        for contour in self.contours3D:
            for primitive in contour:
                primitive.MPLPlot(ax)
        
    def FreeCADExport(self, ip, ndigits=3):
        name = 'primitive'+str(ip)
        s = 'W=[]\n'
        for ic, contour in enumerate(self.contours3d): 
            s += 'L=[]\n'
            for ip, primitive in enumerate(contour.basis_primitives):
                s += primitive.FreeCADExport('L{}_{}'.format(ic,ip),8)
                s += 'L.append(L{}_{})\n'.format(ic,ip)
            s += 'S = Part.Shape(L)\n' 
            s += 'W.append(Part.Wire(S.Edges))\n'
        s += 'F=Part.Face(W)\n'
        a1, a2, a3 = self.axis.vector
        ap1, ap2, ap3 = self.axis_point.vector
        ap1 = round(ap1*1000,ndigits)
        ap2 = round(ap2*1000,ndigits)
        ap3 = round(ap3*1000,ndigits)
        angle = self.angle/math.pi*180
        s += '{} = F.revolve(fc.Vector({},{},{}),fc.Vector({},{},{}),{})\n'.format(name,ap1,ap2,ap3,a1,a2,a3,angle)

#            myObject.Shape = Sweep
        return s
    
    def Volume(self):
        areas=[c.Area() for c in self.contours2D]
        # Maximum area is main surface, others cut into it
        sic=list(npy.argsort(areas))[::-1]# sorted indices of contours
        p1=self.axis_point.PlaneProjection(self.plane_origin,self.x,self.y)
        if self.axis_point.PointDistance(p1)!=0:
            raise NotImplementedError
        p1_2D=p1.To2D(self.axis_point,self.x,self.y)
        p2_3D=self.axis_point+volmdlr.Point3D(self.axis.vector)
        p2=p2_3D.PlaneProjection(self.plane_origin,self.x,self.y)
        if p2_3D.PointDistance(p2)!=0:
            raise NotImplementedError
        p2_2D=p2_3D.To2D(self.plane_origin,self.x,self.y)
        axis_2D=volmdlr.Line2D(p1_2D,p2_2D)
        
        com=self.contours2D[sic[0]].CenterOfMass()
        rg=axis_2D.PointDistance(com)
        volume=areas[sic[0]]*rg
        
        for i in sic[1:]:
            com=self.contours2D[i].CenterOfMass()
            rg=axis_2D.PointDistance(com)
            volume-=areas[i]*rg
            
        return self.angle*volume
                
        
        
class HelicalExtrudedProfile(volmdlr.Primitive3D):
    """

    """
    def __init__(self, plane_origin, x, y, axis_point,axis, pitch,
                 outer_contour2d, inner_contours2d=[], name=''):
        volmdlr.Primitive3D.__init__(self, name)
        self.inner_contours2d = inner_contours2d
        self.outer_contour2d = outer_contour2d
        self.axis_point=axis_point
        self.axis=axis
        self.pitch=pitch
       
        self.inner_contours3d=[c.To3D(plane_origin,x,y) for c in inner_contours2d]
        self.outer_contour3d=outer_contour2d.To3D(plane_origin,x,y)
            
            
    def FreeCADExport(self,ip,ndigits=3):
        name='primitive{}'.format(ip)
        s="E=[]\n"   
        for ip,primitive in enumerate(self.outer_contour3d.basis_primitives):
            s+=primitive.FreeCADExport('L_{}'.format(ip))
            s+='E.append(Part.Edge(L_{}))\n'.format(ip)
#        s+='S = Part.Shape(L)\n' 
        s+='W=Part.Wire(E[:])\n'

        a1,a2,a3=self.axis
        ap1,ap2,ap3=self.axis_point
        ap1=round(ap1*1000,ndigits)
        ap2=round(ap2*1000,ndigits)
        ap3=round(ap3*1000,ndigits)
        
        width = self.axis.Norm()*1000
        direction = bool(self.pitch < 0)
        pitch = round(abs(self.pitch)*1000, ndigits)
        s += "helix = Part.makeHelix({}, {}, 50., 0, {})\n".format(pitch, width, direction)
        s += "helix.translate(fc.Vector({},{},{}))\n".format(ap1, ap2, ap3)

        s += '{} = helix.makePipeShell([W],True,True)\n'.format(name)
        for ic, contour in enumerate(self.inner_contours3d):
            s += "Ei=[]\n"   
            s += "helix2 = Part.makeHelix({}, {}, 50., 0, {})\n".format(pitch,1.01*width+pitch,direction)
            s += "helix2.translate(fc.Vector({},{},{}))\n".format(ap1,ap2,ap3-pitch)
            for ip, primitive in enumerate(contour.basis_primitives):
                s += primitive.FreeCADExport('L_{}_{}'.format(ic,ip))
                s += 'Ei.append(Part.Edge(L_{}_{}))\n'.format(ic,ip)
            s+= 'Wi = Part.Wire(Ei[:])\n'
            s+= "{} = {}.cut(helix2.makePipeShell([Wi],True,True))\n".format(name,name)

        return s
    
    
class Sweep(volmdlr.Primitive3D):
    """
    Sweep a 3D contour along a Wire3D
    """
    def __init__(self, contour3d, wire3d, name=''):
        volmdlr.Primitive3D.__init__(self,name)
        self.contour3d = contour3d
        self.wire3d = wire3d
            
    def FreeCADExport(self,ip,ndigits=3):
        name='primitive{}'.format(ip)
        s="E=[]\n"   
        for ip, primitive in enumerate(self.contour3d.basis_primitives):
            s+=primitive.FreeCADExport('L_{}'.format(ip))
            s+='E.append(Part.Edge(L_{}))\n'.format(ip)
#        s+='S = Part.Shape(L)\n' 
        s+='contour=Part.Wire(E[:])\n'

        s+="E=[]\n"   
        for ip, primitive in enumerate(self.wire3d.basis_primitives):
            s+=primitive.FreeCADExport('L_{}'.format(ip))
            s+='E.append(Part.Edge(L_{}))\n'.format(ip)
#        s+='S = Part.Shape(L)\n' 
        s+='wire=Part.Wire(E[:])\n'

        s+='{}=wire.makePipeShell([contour],True, True)\n'.format(name)


        return s
    
class Cut(volmdlr.Primitive3D):
    """
    Cut primitive 1 by primitive 2
    """
    def __init__(self,primitive,cut_primitive,name=''):
        volmdlr.Primitive3D.__init__(self,name)
        self.primitive=primitive
        self.cut_primitive=cut_primitive
        
    # TODO: this one is oubviously false
#    def Volume(self):
#        return self.primitive.Volume() - self.cut_primitive.Volume()

        
    def FreeCADExport(self,ip):
        name = 'primitive{}'.format(ip)
        
        s = self.primitive.FreeCADExport('{}_0'.format(ip))
        s += self.cut_primitive.FreeCADExport('{}_1'.format(ip))
        
        s+="{} = {}_0.cut({}_1)\n".format(name,name,name)

        return s
        