#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 4 2022

@author: s.bendjebla
"""

import volmdlr 
import volmdlr.core 

# %% VolumeModel: dict_to_object

dict_obj = {
  "name": "",
  "object_class": "volmdlr.core.VolumeModel",
  "package_version": "0.4.1.dev487+g91fdcd3a",
  "primitives": [
    {
      "name": "bounding_box",
      "object_class": "volmdlr.primitives3d.Block",
      "package_version": "0.4.1.dev487+g91fdcd3a",
      "frame": {
        "object_class": "volmdlr.Frame3D",
        "name": "",
        "origin": {
          "object_class": "volmdlr.Point3D",
          "x": 0.19688422265164585,
          "y": 0.0,
          "z": 0.3142481511462552,
          "name": ""
        },
        "u": {
          "object_class": "volmdlr.Vector3D",
          "x": 0.6937684453032917,
          "y": 0.0,
          "z": 0.0,
          "name": ""
        },
        "v": {
          "object_class": "volmdlr.Vector3D",
          "x": 0.0,
          "y": 0.2552053569512336,
          "z": 0.0,
          "name": ""
        },
        "w": {
          "object_class": "volmdlr.Vector3D",
          "x": 0.0,
          "y": 0.0,
          "z": 0.593794544922806,
          "name": ""
        }
      },
      "alpha": 0.9,
      "color": [
        0.8,
        0.8,
        0.8
      ]
    },
    {
      "name": "sensor_free_0",
      "object_class": "volmdlr.primitives3d.Block",
      "package_version": "0.4.1.dev487+g91fdcd3a",
      "frame": {
        "object_class": "volmdlr.Frame3D",
        "name": "",
        "origin": {
          "object_class": "volmdlr.Point3D",
          "x": 0.5437684453032917,
          "y": 0.08361153520018745,
          "z": 0.4786331900941267,
          "name": ""
        },
        "u": {
          "object_class": "volmdlr.Vector3D",
          "x": 0.2,
          "y": 0.0,
          "z": 0.0,
          "name": ""
        },
        "v": {
          "object_class": "volmdlr.Vector3D",
          "x": 0.0,
          "y": 0.2,
          "z": 0.0,
          "name": ""
        },
        "w": {
          "object_class": "volmdlr.Vector3D",
          "x": 0.0,
          "y": 0.0,
          "z": 0.2,
          "name": ""
        }
      },
      "alpha": 1,
      "color": [
        0.2,
        1,
        0.4
      ]
    }
  ]
}

model = volmdlr.core.VolumeModel.dict_to_object(dict_obj)
# m.babylonjs()

# %% geo file generation

model.to_geo('model_two_closed_shell')

# %% gmsh file generation

# model.to_msh('model_two_closed_shell', 2)
