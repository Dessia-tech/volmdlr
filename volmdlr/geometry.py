#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 11:14:32 2017

@author: steven
"""
from numpy import dot,cross,array,zeros,random
import volmdlr as vm
import math
from scipy.linalg import norm


def PointProjectionPlane(point,plane_origin,plane_normal):
    """
    Plane is defined by (plane_point,plane_normal)
    :returns : coordinates in global coordinate system
    """
    return vm.Point2D(point-dot(point-plane_origin,plane_normal)*plane_normal)

def PointLocalProjectionPlane(point,plane_origin,x1,x2):
    """
    Plane is defined by (plane_point,x1,x2)
    :returns : coordinates in local coordinate system
    """
    x1=x1.vector
    x2=x2.vector
    point=point.vector
    plane_origin=plane_origin.vector
    plane_normal=cross(x1,x2)
    xp=point-dot(point-plane_origin,plane_normal)*plane_normal-plane_origin# projeted point
    return vm.Point2D((dot(xp,x1),dot(xp,x2)))

def Euler2TransferMatrix(psi,theta,phi):
    """
    Give Transition Matrix from euler angles
    Angles in radians
    """

    cpsi=math.cos(psi)
    spsi=math.sin(psi)
    ctheta=math.cos(theta)
    stheta=math.sin(theta)
    cphi=math.cos(phi)
    sphi=math.sin(phi)
    P=array([[cphi*cpsi-sphi*ctheta*spsi,-spsi*cphi-cpsi*ctheta*sphi,stheta*sphi],[cpsi*sphi+spsi*ctheta*cphi,-sphi*spsi+cphi*ctheta*cpsi,-stheta*cphi],[spsi*stheta,cpsi*stheta,ctheta]])
    return P
    

def TransferMatrix2Euler(R):
    if ((R[2,2]!=1)&(R[2,2]!=-1)):
        theta=math.acos(R[2,2])
        psi=math.atan2(R[2,0]/math.sin(theta),R[2,1]/math.sin(theta))
        phi=math.atan2(R[0,2]/math.sin(theta),-R[1,2]/math.sin(theta))
    else:
        phi=0
        if R[2,2]==1:
            theta=0
            psi=math.atan2(R[1,0],R[0,0])
        else:
            theta=math.pi
            psi=-math.atan2(R[1,0],R[0,0])
    return(array([psi,theta,phi]))

    
def Direction2Euler(u,v=random.random(3)):
#    u=npy.array([ux,uy,uz])
    u=u/norm(u)
    R=zeros((3,3))
    R[:,0]=u
    v=v-dot(u,v)*u
    v=v/norm(v)
    w=cross(u,v)
    R[:,1]=v
    R[:,2]=w
    euler=TransferMatrix2Euler(R)
    return euler


def Huygens2D(I1,area,point1,point2):
    """
    area acts the same way as the mass in 3D
    """
    a,b=point1.vector-point2.vector
    I2=I1+area*array([[b**2,-a*b],[-a*b,a**2]])
    return I2

