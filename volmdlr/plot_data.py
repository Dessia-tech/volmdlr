#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 14:07:37 2017

@author: steven
"""

import math
import numpy as npy
npy.seterr(divide='raise')
import volmdlr as vm
#from itertools import permutations
import jsonschema

import matplotlib.pyplot as plt
from matplotlib.patches import Arc, FancyArrow
from mpl_toolkits.mplot3d import Axes3D

_jsonschema = {
    "definitions": {},
    "$schema": "http://json-schema.org/draft-07/schema#",
    "title": "plot_data Base Schema",
    'type' : 'array',
    'editable' : True,
    'description' : 'plot data',
    "items" : {
          "anyOf": [
               {"type" : "object",
                "required": ['type', 'data', 'color', 'marker', 'size', 'dash'],
                "properties": {'type' : {"const": 'line'},
                               'data' : {'type' : 'array',
                                         "minItems": 4,
                                         "maxItems": 4,
                                         'examples' : [[0, 1, 1, 2]],
                                         "description" : "Point coord",
                                         'items' : {'type' : 'number',
                                                    'step' : "any",
                                                    'unit' : 'm'
                                                    }
                                         },
                                'color': {"enum": ['black', 'red', 'blue']},
                                'marker': {"enum": [None, 'o']},
                                'dash': {'type': 'boolean'},
                                'size': {"type": 'number',
                                         "examples" : [0.5],
                                         "step" : 'any',
                                         "minimum" : 0,
                                         "description" : "Point size"},
                                'opacity': {"type": 'number',
                                            "examples" : [1],
                                            "step" : 'any',
                                            "minimum" : 0,
                                            "description" : "Point opacity"},
                               }
                       },
               {"type" : "object",
                "required": ['type', 'data', 'color', 'marker', 'size'],
                "properties": {'type' : {"const": 'point'},
                               'data' : {'type' : 'array',
                                         "minItems": 2,
                                         "maxItems": 2,
                                         'examples' : [[0, 1]],
                                         "description" : "Point coord",
                                         'items' : {'type' : 'number',
                                                    'step' : "any",
                                                    'unit' : 'm'
                                                    }
                                         },
                                'color': {"enum": ['black', 'red', 'blue']},
                                'marker': {"enum": ['o']},
                                'size': {"type": 'number',
                                         "examples" : [0.5],
                                         "step" : 'any',
                                         "minimum" : 0,
                                         "description" : "Point size"},
                                'opacity': {"type": 'number',
                                            "examples" : [1],
                                            "step" : 'any',
                                            "minimum" : 0,
                                            "description" : "Point opacity"},
                               }
                       },
               {"type" : "object",
                "properties": {'type' : {"const": 'contour'}}
                       },
               {"type" : "object",
                "properties": {'type' : {"const": 'wire'}}
                       },
                    ]
                    },
         }
     
         
color = {'black': 'k', 'blue': 'b', 'red': 'r'}

def validate(plot_datas):
    return jsonschema.validate(instance=plot_datas, schema=_jsonschema)

def plot(plot_datas, ax=None):
    if ax is None:
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
    else:
        fig = None
        
    for plot_data in plot_datas:
        if plot_data['type'] == 'line':
            style = ''
            if plot_data['dash']:
                style += '--'
            else:
                style += '-'
            style += color[plot_data['color']]
            p1, p2 = plot_data['data'][0: 2], plot_data['data'][2:] 
            if plot_data['arrow']:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], style)
                length = ((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)**0.5
                if width is None:
                    width = length / 1000.
                    head_length = length/20.
                    head_width = head_length/2.
                else:
                    head_width = 2*width
                    head_length = head_width
                ax.arrow(p1[0], p1[1], (p2[0] - p1[0])/length*(length - head_length),
                         (p2[1] - p1[1])/length*(length - head_length),
                         head_width = head_width, fc = 'b', linewidth = 0,
                         head_length = head_length, width = width, alpha = 0.3)
            else:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], style, linewidth=plot_data['size'])
                
        elif plot_data['type'] == 'point':
            p1 = plot_data['data']
            style = ''
            style += color[plot_data['color']]
            style += plot_data['marker']
            ax.plot(p1[0], p1[1], style, linewidth=plot_data['size'])
            
        elif plot_data['type'] == 'contour':
            plot(plot_data['plot_data'], ax)
            
        elif plot_data['type'] == 'arc':
            pc = vm.Point2D(( plot_data['cx'],  plot_data['cy']))
            ax.add_patch(Arc(pc, 2*plot_data['r'], 2*plot_data['r'], angle=0,
                        theta1=plot_data['angle1']*0.5/math.pi*360,
                        theta2=plot_data['angle2']*0.5/math.pi*360,
                        color=color[plot_data['color']], linewidth=plot_data['size']))
            
        elif plot_data['type'] == 'circle':
            pc = vm.Point2D(( plot_data['cx'],  plot_data['cy']))
            ax.add_patch(Arc(pc, 2*plot_data['r'], 2*plot_data['r'], angle=0,
                        theta1=0,
                        theta2=360,
                        color=color[plot_data['color']], linewidth=plot_data['size']))
    return fig, ax

    