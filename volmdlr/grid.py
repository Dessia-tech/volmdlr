#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Module containing grid and relative objects
"""

import volmdlr as vm
import volmdlr.wires
import numpy as npy
from dessia_common import DessiaObject



class Grid2D(DessiaObject):
    
    def __init__(self, list_points: list[list[volmdlr.Point2D]],
                 direction: list[str] = ['+x','+y'],
                 name: str = ''):
        self.list_points = list_points
        direction = direction
        DessiaObject.__init__(self, name=name)
