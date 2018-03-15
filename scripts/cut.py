#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 12:34:01 2018

@author: Steven Masfaraud masfaraud@dessia.tech
"""

import volmdlr as vm
import volmdlr.primitives3D as p3D

sphere=p3D.Sphere(vm.Point3D((0,0,0)),4)

cylinder=p3D.Cylinder((0,0,0),(0,0,1),0.5,10)

cut=p3D.Cut(sphere,cylinder)

model=vm.VolumeModel([cut])

model.FreeCADExport('python','cut','/usr/lib/freecad/lib')