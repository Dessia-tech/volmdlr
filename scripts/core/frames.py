#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 15:35:13 2022

@author: steven
"""

import volmdlr as vm

frame1 = vm.OXYZ
center = vm.Point3D(0.6, 0, 0)
frame2 = frame1.rotation(center, vm.Z3D, 0.2)

ax = frame1.plot()
center.plot(ax=ax, color='r')
frame2.plot(ax=ax)