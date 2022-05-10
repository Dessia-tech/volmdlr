#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 17:02:36 2017

"""

import math
import volmdlr as vm
import volmdlr.wires as vmw
import time



p1 = vm.Point3D(0, 0, 0)
p2 = vm.Point3D(1, 0, 0.1)
p3 = vm.Point3D(2, 1, -0.1)
p4 = vm.Point3D(1, 0.5, 0.1)
p5 = vm.Point3D(-0.5, 1, 0)
polygon = vm.wires.ClosedPolygon3D([p1, p2, p3, p4, p5])


polygon._check_platform()