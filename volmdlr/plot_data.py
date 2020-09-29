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
# from itertools import permutations
import jsonschema
import json

import matplotlib.pyplot as plt
from matplotlib.patches import Arc, FancyArrow
from mpl_toolkits.mplot3d import Axes3D
import pkg_resources
import tempfile
import webbrowser
from dessia_common import DessiaObject
from typing import TypeVar, List

from jinja2 import Environment, PackageLoader, select_autoescape


class ColorMapSet(DessiaObject):
    def __init__(self, value: float = None, tooltip: bool = False,
                 color_range: str = None, selector: bool = True,
                 name: str = ''):
        self.selector = selector
        self.color_range = color_range
        self.tooltip = tooltip
        self.value = value
        DessiaObject.__init__(self, name=name)


class HatchingSet(DessiaObject):
    def __init__(self, stroke_width: float = 1, hatch_spacing: float = 10,
                 name: str = ''):
        self.stroke_width = stroke_width
        self.hatch_spacing = hatch_spacing
        DessiaObject.__init__(self, name=name)


class ColorSurfaceSet(DessiaObject):
    def __init__(self, color: str = 'white', name: str = ''):
        self.color = color
        DessiaObject.__init__(self, name=name)


class PointShapeSet(DessiaObject):
    def __init__(self, shape: str = 'circle', name: str = ''):
        self.shape = shape
        DessiaObject.__init__(self, name=name)


class PointSizeSet(DessiaObject):
    def __init__(self, size: int, name: str = ''):
        self.size = size
        DessiaObject.__init__(self, name=name)


class PointColorSet(DessiaObject):
    def __init__(self, color_fill: str, color_stroke: str, name: str = ''):
        self.color_fill = color_fill
        self.color_stroke = color_stroke
        DessiaObject.__init__(self, name=name)


class WindowSizeSet(DessiaObject):
    def __init__(self, width: float, height: float, name: str = ''):
        self.width = width
        self.height = height
        DessiaObject.__init__(self, name=name)


class PlotDataState(DessiaObject):
    def __init__(self, name: str = '', color_map: ColorMapSet = None,
                 hatching: HatchingSet = None,
                 color_surface: ColorSurfaceSet = None,
                 shape_set: PointShapeSet = None,
                 point_size: PointSizeSet = None,
                 point_color: PointColorSet = None,
                 window_size: WindowSizeSet = None,
                 stroke_width: float = 1, color_line: str = 'black',
                 marker: str = None,
                 dash: str = None, opacity: float = 1):
        self.color_surface = color_surface
        self.color_map = color_map
        self.hatching = hatching
        self.opacity = opacity
        self.dash = dash
        self.marker = marker
        self.color_line = color_line
        self.stroke_width = stroke_width
        self.shape_set = shape_set
        if self.shape_set is None:
            self.shape_set = PointShapeSet(shape='circle')
        self.point_size = point_size
        if self.point_size is None:
            self.point_size = PointSizeSet(size=2)
        self.point_color = point_color
        if self.point_color is None:
            self.point_color = PointColorSet(color_fill='black',
                                             color_stroke='black')
        self.window_size = window_size
        DessiaObject.__init__(self, name=name)


class PlotDataLine2D(DessiaObject):
    def __init__(self, data: List[float],
                 plot_data_states: List[PlotDataState],
                 type: str = 'line', name: str = '', ):
        self.data = data
        self.type = type
        self.plot_data_states = plot_data_states
        DessiaObject.__init__(self, name=name)


class PlotDataCircle2D(DessiaObject):
    def __init__(self, cx: float, cy: float, r: float,
                 plot_data_states: List[PlotDataState],
                 type: str = 'circle', name: str = '', ):
        self.type = type
        self.plot_data_states = plot_data_states
        self.r = r
        self.cy = cy
        self.cx = cx
        DessiaObject.__init__(self, name=name)


class PlotDataPoint2D(DessiaObject):
    def __init__(self, cx: float, cy: float,
                 plot_data_states: List[PlotDataState], type: str = 'point',
                 name: str = '', ):
        self.type = type
        self.plot_data_states = plot_data_states
        self.cx = cx
        self.cy = cy
        DessiaObject.__init__(self, name=name)


class PlotDataScatterPlot(DessiaObject):
    def __init__(self, nb_points_x: int, nb_points_y: int, font_size: int,
                 graduation_color: str, axis_color: str, plot_data_states: List[PlotDataState],
                 arrow_on: bool, axis_width:str, grid_on: bool,
                 name: str = '',
                 type: str = 'axis'):
        self.nb_points_x = nb_points_x
        self.nb_points_y = nb_points_y
        self.font_size = font_size
        self.graduation_color = graduation_color
        self.axis_color = axis_color
        self.plot_data_states = plot_data_states
        self.arrow_on = arrow_on
        self.axis_width = axis_width
        self.grid_on = grid_on
        self.type = type
        DessiaObject.__init__(self, name=name)


class PlotDataTooltip(DessiaObject):
    def __init__(self, colorfill: str, font: str, tp_width: float, tp_radius: float, to_plot_list: list,
                 plot_data_states: List[PlotDataState], type: str = 'tooltip',
                 name: str = ''):
        self.colorfill = colorfill
        self.font = font
        self.tp_width = tp_width
        self.tp_radius = tp_radius
        self.to_plot_list = to_plot_list
        self.plot_data_states = plot_data_states
        self.type = type
        DessiaObject.__init__(self, name=name)


class PlotDataArc2D(DessiaObject):
    def __init__(self, cx: float, cy: float, r: float,
                 data: List[float], angle1: float, angle2: float,
                 plot_data_states: List[PlotDataState],
                 type: str = 'arc', name: str = '', ):
        self.angle2 = angle2
        self.angle1 = angle1
        self.data = data
        self.type = type
        self.plot_data_states = plot_data_states
        self.r = r
        self.cy = cy
        self.cx = cx
        DessiaObject.__init__(self, name=name)


class PlotDataContour2D(DessiaObject):
    def __init__(self, plot_data_primitives: List[float],
                 plot_data_states: List[PlotDataState],
                 type: str = 'contour', name: str = '', ):
        self.plot_data_primitives = plot_data_primitives
        self.type = type
        self.plot_data_states = plot_data_states
        DessiaObject.__init__(self, name=name)


color = {'black': 'k', 'blue': 'b', 'red': 'r', 'green': 'g'}


def plot_d3(plot_datas):
    env = Environment(loader=PackageLoader('volmdlr', 'templates'),
                      autoescape=select_autoescape(['html', 'xml']))
    template = env.get_template('plot_data3.html')

    volmdlr_path = pkg_resources.resource_filename(
        pkg_resources.Requirement('volmdlr'),
        'volmdlr/templates')
    data = []
    for d in plot_datas:
        data.append(json.dumps(d))
    s = template.render(
        volmdlr_path=volmdlr_path,
        D3Data=data)
    temp_file = tempfile.mkstemp(suffix='.html')[1]

    with open(temp_file, 'wb') as file:
        file.write(s.encode('utf-8'))

    webbrowser.open('file://' + temp_file)
    print('file://' + temp_file)


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
                length = ((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2) ** 0.5
                if width is None:
                    width = length / 1000.
                    head_length = length / 20.
                    head_width = head_length / 2.
                else:
                    head_width = 2 * width
                    head_length = head_width
                ax.arrow(p1[0], p1[1],
                         (p2[0] - p1[0]) / length * (length - head_length),
                         (p2[1] - p1[1]) / length * (length - head_length),
                         head_width=head_width, fc='b', linewidth=0,
                         head_length=head_length, width=width, alpha=0.3)
            else:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], style,
                        linewidth=plot_data['size'])

        elif plot_data['type'] == 'point':
            p1 = plot_data['data']
            style = ''
            style += color[plot_data['color']]
            style += plot_data['marker']
            ax.plot(p1[0], p1[1], style, linewidth=plot_data['size'])

        elif plot_data['type'] == 'contour':
            plot(plot_data['plot_data'], ax)

        elif plot_data['type'] == 'arc':
            pc = vm.Point2D((plot_data['cx'], plot_data['cy']))
            ax.add_patch(
                Arc(pc, 2 * plot_data['r'], 2 * plot_data['r'], angle=0,
                    theta1=plot_data['angle1'] * 0.5 / math.pi * 360,
                    theta2=plot_data['angle2'] * 0.5 / math.pi * 360,
                    color=color[plot_data['color']],
                    linewidth=plot_data['size']))

        elif plot_data['type'] == 'circle':
            pc = vm.Point2D((plot_data['cx'], plot_data['cy']))
            ax.add_patch(
                Arc(pc, 2 * plot_data['r'], 2 * plot_data['r'], angle=0,
                    theta1=0,
                    theta2=360,
                    color=color[plot_data['color']],
                    linewidth=plot_data['size']))
    return fig, ax
