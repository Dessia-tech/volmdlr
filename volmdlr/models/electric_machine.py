#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 17:26:59 2020

@author: mgasmi
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import volmdlr.edges as edges
import volmdlr.wires as wires
import volmdlr.faces as faces
import volmdlr.core_compiled as vm
import volmdlr.mesh as vmmesh

rotor_external=wires.Contour2D.dict_to_object({"name": "","object_class": "volmdlr.wires.Contour2D","package_version": "0.1.13.dev167","primitives": [{"name": "","object_class": "volmdlr.edges.LineSegment2D","package_version": "0.1.13.dev167","points": [{"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.05,"y": 0.055999999999999994},{"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.05,"y": 0.055999999999999994}],"start": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.05,"y": 0.055999999999999994},"end": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.05,"y": 0.055999999999999994}},{"name": "","object_class": "volmdlr.edges.LineSegment2D","package_version": "0.1.13.dev167","points": [{"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.05,"y": 0.055999999999999994},{"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.055999999999999994,"y": 0.05}],"start": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.05,"y": 0.055999999999999994},"end": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.055999999999999994,"y": 0.05}},{"name": "","object_class": "volmdlr.edges.LineSegment2D","package_version": "0.1.13.dev167","points": [{"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.055999999999999994,"y": 0.05},{"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.055999999999999994,"y": -0.05}],"start": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.055999999999999994,"y": 0.05},"end": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.055999999999999994,"y": -0.05}},{"name": "","object_class": "volmdlr.edges.LineSegment2D","package_version": "0.1.13.dev167","points": [{"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.055999999999999994,"y": -0.05},{"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.05000000000000001,"y": -0.05599999999999999}],"start": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.055999999999999994,"y": -0.05},"end": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.05000000000000001,"y": -0.05599999999999999}},{"name": "","object_class": "volmdlr.edges.LineSegment2D","package_version": "0.1.13.dev167","points": [{"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.05000000000000001,"y": -0.05599999999999999},{"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.049999999999999996,"y": -0.056}],"start": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.05000000000000001,"y": -0.05599999999999999},"end": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.049999999999999996,"y": -0.056}},{"name": "","object_class": "volmdlr.edges.LineSegment2D","package_version": "0.1.13.dev167","points": [{"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.049999999999999996,"y": -0.056},{"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.05599999999999999,"y": -0.05000000000000001}],"start": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.049999999999999996,"y": -0.056},"end": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.05599999999999999,"y": -0.05000000000000001}},{"name": "","object_class": "volmdlr.edges.LineSegment2D","package_version": "0.1.13.dev167","points": [{"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.05599999999999999,"y": -0.05000000000000001},{"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.056,"y": 0.049999999999999996}],"start": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.05599999999999999,"y": -0.05000000000000001},"end": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.056,"y": 0.049999999999999996}},{"name": "","object_class": "volmdlr.edges.LineSegment2D","package_version": "0.1.13.dev167","points": [{"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.056,"y": 0.049999999999999996},{"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.05,"y": 0.055999999999999994}],"start": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.056,"y": 0.049999999999999996},"end": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.05,"y": 0.055999999999999994}}]})

rotor_internal=wires.Contour2D.dict_to_object({"name": "","object_class": "volmdlr.wires.Contour2D","package_version": "0.1.13.dev167","primitives": [{"name": "","object_class": "volmdlr.wires.Circle2D","package_version": "0.1.13.dev167","center": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.0,"y": 0.0},"radius": 0.03,"angle": 6.283185307179586}]})

rotor_magnet=wires.Contour2D.dict_to_object({"name": "","object_class": "volmdlr.wires.Contour2D","package_version": "0.1.13.dev167","primitives": [{"name": "","object_class": "volmdlr.edges.Arc2D","package_version": "0.1.13.dev167","start": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.05,"y": 0.08077747210701756},"end": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.05,"y": 0.08077747210701756},"interior": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.0,"y": 0.095},"center": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 4.336808689942018e-18,"y": -4.573879611178527e-17},"radius": 0.09500000000000004,"is_trigo": False,"angle1": 1.0165344923425688,"angle2": 2.1250581612472246,"angle": 1.1085236689046558},{"name": "","object_class": "volmdlr.edges.LineSegment2D","package_version": "0.1.13.dev167","points": [{"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.05,"y": 0.08077747210701756},{"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.05,"y": 0.055999999999999994}],"start": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.05,"y": 0.08077747210701756},"end": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.05,"y": 0.055999999999999994}},{"name": "","object_class": "volmdlr.edges.LineSegment2D","package_version": "0.1.13.dev167","points": [{"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.05,"y": 0.055999999999999994},{"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.05,"y": 0.055999999999999994}],"start": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.05,"y": 0.055999999999999994},"end": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.05,"y": 0.055999999999999994}},{"name": "","object_class": "volmdlr.edges.LineSegment2D","package_version": "0.1.13.dev167","points": [{"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.05,"y": 0.055999999999999994},{"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.05,"y": 0.08077747210701756}],"start": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.05,"y": 0.055999999999999994},"end": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": -0.05,"y": 0.08077747210701756}}]})                         

stator_external=wires.Contour2D.dict_to_object({"name": "","object_class": "volmdlr.wires.Contour2D","package_version": "0.1.13.dev167","primitives": [{"name": "","object_class": "volmdlr.wires.Circle2D","package_version": "0.1.13.dev167","center": {"name": "","object_class": "volmdlr.core_compiled.Point2D","package_version": "0.1.13.dev167","x": 0.0,"y": 0.0},"radius": 0.2,"angle": 6.283185307179586}]})

stator_internal=wires.Contour2D.dict_to_object({
"name": "",
"object_class": "volmdlr.wires.Contour2D",
"package_version": "0.1.13.dev167",
"primitives": [
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.01,
"y": 0.1
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.01,
"y": 0.115
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.01,
"y": 0.1
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.01,
"y": 0.115
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.01,
"y": 0.115
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.02176493326044296,
"y": 0.11793301815030346
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.01,
"y": 0.115
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.02176493326044296,
"y": 0.11793301815030346
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.02176493326044296,
"y": 0.11793301815030346
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.035791046509173426,
"y": 0.17027918542728054
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.02176493326044296,
"y": 0.11793301815030346
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.035791046509173426,
"y": 0.17027918542728054
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.035791046509173426,
"y": 0.17027918542728054
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.035791046509173426,
"y": 0.17027918542728054
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.0,
"y": 0.17400000000000002
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.0,
"y": 0.0
},
"radius": 0.17400000000000002,
"is_trigo": False,
"angle1": 1.3636218017252826,
"angle2": 1.7779708518645108,
"angle": 0.4143490501392282
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.035791046509173426,
"y": 0.17027918542728054
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.02176493326044296,
"y": 0.11793301815030346
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.035791046509173426,
"y": 0.17027918542728054
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.02176493326044296,
"y": 0.11793301815030346
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.02176493326044296,
"y": 0.11793301815030346
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.01,
"y": 0.115
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.02176493326044296,
"y": 0.11793301815030346
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.01,
"y": 0.115
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.01,
"y": 0.115
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.01,
"y": 0.1
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.01,
"y": 0.115
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.01,
"y": 0.1
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.01,
"y": 0.1
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.04133974596215561,
"y": 0.09160254037844388
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.01622264624736139,
"y": 0.09918077307993205
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -2.6733053435393526e-17,
"y": -2.030577005775927e-16
},
"radius": 0.10049875621120913,
"is_trigo": False,
"angle1": 1.1468662036877606,
"angle2": 1.4711276743037345,
"angle": 0.3242614706159739
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.04133974596215561,
"y": 0.09160254037844388
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.04883974596215561,
"y": 0.10459292143521046
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.04133974596215561,
"y": 0.09160254037844388
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.04883974596215561,
"y": 0.10459292143521046
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.04883974596215561,
"y": 0.10459292143521046
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.04011752395993525,
"y": 0.11301545629335558
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.04883974596215561,
"y": 0.10459292143521046
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.04011752395993525,
"y": 0.11301545629335558
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.04011752395993525,
"y": 0.11301545629335558
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.05414363720866572,
"y": 0.16536162357033266
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.04011752395993525,
"y": 0.11301545629335558
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.05414363720866572,
"y": 0.16536162357033266
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.05414363720866572,
"y": 0.16536162357033266
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.1161355482186148,
"y": 0.12957057706115924
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.087,
"y": 0.15068842025849236
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 3.2867877316145964e-16,
"y": 6.177564192484373e-16
},
"radius": 0.17399999999999932,
"is_trigo": False,
"angle1": 0.8400230261269829,
"angle2": 1.2543720762662127,
"angle": 0.41434905013922985
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.1161355482186148,
"y": 0.12957057706115924
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.0778154941903682,
"y": 0.09125052303291262
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.1161355482186148,
"y": 0.12957057706115924
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.0778154941903682,
"y": 0.09125052303291262
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.0778154941903682,
"y": 0.09125052303291262
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.06616025403784438,
"y": 0.09459292143521045
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.0778154941903682,
"y": 0.09125052303291262
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.06616025403784438,
"y": 0.09459292143521045
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.06616025403784438,
"y": 0.09459292143521045
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.05866025403784438,
"y": 0.08160254037844387
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.06616025403784438,
"y": 0.09459292143521045
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.05866025403784438,
"y": 0.08160254037844387
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.05866025403784438,
"y": 0.08160254037844387
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.08160254037844386,
"y": 0.05866025403784441
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.06363961030678927,
"y": 0.07778174593052024
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -1.2468035368457382e-16,
"y": -1.6248660041655238e-16
},
"radius": 0.10049875621120911,
"is_trigo": False,
"angle1": 0.6232674280894617,
"angle2": 0.9475288987054358,
"angle": 0.3242614706159741
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.08160254037844386,
"y": 0.0586602540378444
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.09459292143521043,
"y": 0.0661602540378444
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.08160254037844386,
"y": 0.0586602540378444
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.09459292143521043,
"y": 0.0661602540378444
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.09459292143521043,
"y": 0.0661602540378444
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.0912505230329126,
"y": 0.07781549419036822
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.09459292143521043,
"y": 0.0661602540378444
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.0912505230329126,
"y": 0.07781549419036822
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.0912505230329126,
"y": 0.07781549419036822
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.1295705770611592,
"y": 0.11613554821861483
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.0912505230329126,
"y": 0.07781549419036822
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.1295705770611592,
"y": 0.11613554821861483
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.1295705770611592,
"y": 0.11613554821861483
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.16536162357033263,
"y": 0.054143637208665746
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.15068842025849233,
"y": 0.08700000000000002
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 8.276436788999663e-17,
"y": 8.975862275846609e-17
},
"radius": 0.1739999999999999,
"is_trigo": False,
"angle1": 0.31642425052868456,
"angle2": 0.730773300667913,
"angle": 0.41434905013922846
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.16536162357033263,
"y": 0.054143637208665746
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.11301545629335556,
"y": 0.04011752395993527
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.16536162357033263,
"y": 0.054143637208665746
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.11301545629335556,
"y": 0.04011752395993527
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.11301545629335556,
"y": 0.04011752395993527
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.10459292143521044,
"y": 0.04883974596215563
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.11301545629335556,
"y": 0.04011752395993527
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.10459292143521044,
"y": 0.04883974596215563
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.10459292143521044,
"y": 0.04883974596215563
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.09160254037844387,
"y": 0.04133974596215563
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.10459292143521044,
"y": 0.04883974596215563
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.09160254037844387,
"y": 0.04133974596215563
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.09160254037844387,
"y": 0.04133974596215563
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.10000000000000002,
"y": 0.010000000000000023
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.09400439217788162,
"y": 0.03554116277314279
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -8.993521131875374e-16,
"y": -2.790353009353527e-16
},
"radius": 0.10049875621120985,
"is_trigo": False,
"angle1": 0.0996686524911641,
"angle2": 0.4239301231071358,
"angle": 0.32426147061597166
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.1,
"y": 0.010000000000000007
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.115,
"y": 0.010000000000000007
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.1,
"y": 0.010000000000000007
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.115,
"y": 0.010000000000000007
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.115,
"y": 0.010000000000000007
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.11793301815030346,
"y": 0.021764933260442966
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.115,
"y": 0.010000000000000007
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.11793301815030346,
"y": 0.021764933260442966
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.11793301815030346,
"y": 0.021764933260442966
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.17027918542728054,
"y": 0.03579104650917344
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.11793301815030346,
"y": 0.021764933260442966
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.17027918542728054,
"y": 0.03579104650917344
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.17027918542728054,
"y": 0.03579104650917344
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.17027918542728054,
"y": -0.03579104650917341
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.17400000000000002,
"y": 1.0654427152581974e-17
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -1.311247760626477e-16,
"y": 1.3631648839701518e-17
},
"radius": 0.17400000000000015,
"is_trigo": False,
"angle1": -0.20717452506961392,
"angle2": 0.20717452506961392,
"angle": 0.41434905013922785
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.17027918542728054,
"y": -0.03579104650917341
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.11793301815030346,
"y": -0.021764933260442952
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.17027918542728054,
"y": -0.03579104650917341
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.11793301815030346,
"y": -0.021764933260442952
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.11793301815030346,
"y": -0.021764933260442952
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.115,
"y": -0.009999999999999993
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.11793301815030346,
"y": -0.021764933260442952
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.115,
"y": -0.009999999999999993
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.115,
"y": -0.009999999999999993
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.1,
"y": -0.009999999999999993
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.115,
"y": -0.009999999999999993
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.1,
"y": -0.009999999999999993
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.1,
"y": -0.009999999999999993
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.09160254037844388,
"y": -0.0413397459621556
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.09918077307993205,
"y": -0.01622264624736138
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -2.2411551149754653e-16,
"y": 6.761925212014346e-17
},
"radius": 0.10049875621120913,
"is_trigo": False,
"angle1": -0.4239301231071364,
"angle2": -0.0996686524911624,
"angle": 0.324261470615974
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.09160254037844387,
"y": -0.041339745962155595
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.10459292143521046,
"y": -0.04883974596215559
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.09160254037844387,
"y": -0.041339745962155595
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.10459292143521046,
"y": -0.04883974596215559
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.10459292143521046,
"y": -0.04883974596215559
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.11301545629335558,
"y": -0.04011752395993523
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.10459292143521046,
"y": -0.04883974596215559
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.11301545629335558,
"y": -0.04011752395993523
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.11301545629335558,
"y": -0.04011752395993523
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.16536162357033266,
"y": -0.05414363720866569
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.11301545629335558,
"y": -0.04011752395993523
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.16536162357033266,
"y": -0.05414363720866569
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.16536162357033266,
"y": -0.05414363720866569
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.12957057706115924,
"y": -0.11613554821861477
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.15068842025849236,
"y": -0.08699999999999997
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 2.8515786235400805e-16,
"y": -1.4714655497584672e-16
},
"radius": 0.1739999999999997,
"is_trigo": False,
"angle1": -0.7307733006679131,
"angle2": -0.3164242505286842,
"angle": 0.4143490501392289
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.12957057706115924,
"y": -0.11613554821861477
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.09125052303291262,
"y": -0.07781549419036818
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.12957057706115924,
"y": -0.11613554821861477
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.09125052303291262,
"y": -0.07781549419036818
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.09125052303291262,
"y": -0.07781549419036818
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.09459292143521045,
"y": -0.06616025403784437
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.09125052303291262,
"y": -0.07781549419036818
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.09459292143521045,
"y": -0.06616025403784437
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.09459292143521045,
"y": -0.06616025403784437
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.08160254037844389,
"y": -0.05866025403784437
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.09459292143521045,
"y": -0.06616025403784437
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.08160254037844389,
"y": -0.05866025403784437
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.08160254037844389,
"y": -0.05866025403784437
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.05866025403784442,
"y": -0.08160254037844386
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.07778174593052024,
"y": -0.06363961030678926
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -3.0902840638587913e-16,
"y": 2.80673721336895e-16
},
"radius": 0.10049875621120932,
"is_trigo": False,
"angle1": -0.9475288987054346,
"angle2": -0.6232674280894611,
"angle": 0.32426147061597344
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.058660254037844424,
"y": -0.08160254037844385
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.06616025403784442,
"y": -0.09459292143521042
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.058660254037844424,
"y": -0.08160254037844385
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.06616025403784442,
"y": -0.09459292143521042
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.06616025403784442,
"y": -0.09459292143521042
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.07781549419036825,
"y": -0.09125052303291259
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.06616025403784442,
"y": -0.09459292143521042
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.07781549419036825,
"y": -0.09125052303291259
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.07781549419036825,
"y": -0.09125052303291259
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.11613554821861485,
"y": -0.1295705770611592
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.07781549419036825,
"y": -0.09125052303291259
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.11613554821861485,
"y": -0.1295705770611592
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.11613554821861485,
"y": -0.1295705770611592
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.054143637208665794,
"y": -0.16536162357033263
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.08700000000000006,
"y": -0.1506884202584923
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -6.573575463229242e-16,
"y": 1.2355128384968829e-15
},
"radius": 0.1740000000000014,
"is_trigo": False,
"angle1": -1.25437207626621,
"angle2": -0.8400230261269854,
"angle": 0.41434905013922463
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.054143637208665794,
"y": -0.16536162357033263
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.0401175239599353,
"y": -0.11301545629335555
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.054143637208665794,
"y": -0.16536162357033263
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.0401175239599353,
"y": -0.11301545629335555
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.0401175239599353,
"y": -0.11301545629335555
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.04883974596215566,
"y": -0.10459292143521043
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.0401175239599353,
"y": -0.11301545629335555
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.04883974596215566,
"y": -0.10459292143521043
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.04883974596215566,
"y": -0.10459292143521043
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.04133974596215565,
"y": -0.09160254037844386
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.04883974596215566,
"y": -0.10459292143521043
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.04133974596215565,
"y": -0.09160254037844386
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.04133974596215565,
"y": -0.09160254037844386
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.01000000000000005,
"y": -0.10000000000000002
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.03554116277314281,
"y": -0.09400439217788162
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 1.1727162454588165e-16,
"y": -6.442408349218695e-16
},
"radius": 0.10049875621120828,
"is_trigo": False,
"angle1": -1.4711276743037347,
"angle2": -1.146866203687758,
"angle": 0.32426147061597677
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.010000000000000012,
"y": -0.1
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.010000000000000014,
"y": -0.115
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.010000000000000012,
"y": -0.1
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.010000000000000014,
"y": -0.115
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.010000000000000014,
"y": -0.115
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.021764933260442973,
"y": -0.11793301815030346
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.010000000000000014,
"y": -0.115
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.021764933260442973,
"y": -0.11793301815030346
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.021764933260442973,
"y": -0.11793301815030346
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.035791046509173446,
"y": -0.17027918542728054
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.021764933260442973,
"y": -0.11793301815030346
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.035791046509173446,
"y": -0.17027918542728054
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.035791046509173446,
"y": -0.17027918542728054
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.035791046509173405,
"y": -0.17027918542728054
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 2.1308854305163948e-17,
"y": -0.17400000000000002
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 2.4234042381691588e-17,
"y": 2.331107130002626e-16
},
"radius": 0.17400000000000024,
"is_trigo": False,
"angle1": -1.7779708518645105,
"angle2": -1.3636218017252828,
"angle": 0.41434905013922774
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.035791046509173405,
"y": -0.17027918542728054
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.021764933260442945,
"y": -0.11793301815030346
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.035791046509173405,
"y": -0.17027918542728054
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.021764933260442945,
"y": -0.11793301815030346
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.021764933260442945,
"y": -0.11793301815030346
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.009999999999999986,
"y": -0.115
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.021764933260442945,
"y": -0.11793301815030346
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.009999999999999986,
"y": -0.115
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.009999999999999986,
"y": -0.115
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.009999999999999988,
"y": -0.1
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.009999999999999986,
"y": -0.115
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.009999999999999988,
"y": -0.1
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.009999999999999988,
"y": -0.1
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.041339745962155595,
"y": -0.09160254037844388
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.016222646247361375,
"y": -0.09918077307993205
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.0,
"y": 0.0
},
"radius": 0.1004987562112089,
"is_trigo": False,
"angle1": -1.9947264499020332,
"angle2": -1.6704649792860584,
"angle": 0.32426147061597477
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.04133974596215559,
"y": -0.0916025403784439
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.04883974596215558,
"y": -0.10459292143521046
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.04133974596215559,
"y": -0.0916025403784439
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.04883974596215558,
"y": -0.10459292143521046
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.04883974596215558,
"y": -0.10459292143521046
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.04011752395993522,
"y": -0.11301545629335559
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.04883974596215558,
"y": -0.10459292143521046
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.04011752395993522,
"y": -0.11301545629335559
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.04011752395993522,
"y": -0.11301545629335559
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.05414363720866568,
"y": -0.16536162357033268
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.04011752395993522,
"y": -0.11301545629335559
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.05414363720866568,
"y": -0.16536162357033268
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.05414363720866568,
"y": -0.16536162357033268
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.11613554821861477,
"y": -0.12957057706115926
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.08699999999999997,
"y": -0.15068842025849236
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 3.70653365839329e-16,
"y": 5.935223768667494e-16
},
"radius": 0.1740000000000007,
"is_trigo": False,
"angle1": -2.3015696274628086,
"angle2": -1.887220577323582,
"angle": 0.4143490501392266
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.11613554821861477,
"y": -0.12957057706115926
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.07781549419036818,
"y": -0.09125052303291263
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.11613554821861477,
"y": -0.12957057706115926
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.07781549419036818,
"y": -0.09125052303291263
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.07781549419036818,
"y": -0.09125052303291263
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.06616025403784435,
"y": -0.09459292143521048
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.07781549419036818,
"y": -0.09125052303291263
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.06616025403784435,
"y": -0.09459292143521048
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.06616025403784435,
"y": -0.09459292143521048
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.05866025403784436,
"y": -0.08160254037844389
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.06616025403784435,
"y": -0.09459292143521048
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.05866025403784436,
"y": -0.08160254037844389
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.05866025403784436,
"y": -0.08160254037844389
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.08160254037844385,
"y": -0.05866025403784443
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.06363961030678925,
"y": -0.07778174593052026
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 0.0,
"y": 0.0
},
"radius": 0.1004987562112089,
"is_trigo": False,
"angle1": -2.5183252255003317,
"angle2": -2.194063754884357,
"angle": 0.32426147061597455
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.08160254037844385,
"y": -0.05866025403784443
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.09459292143521042,
"y": -0.06616025403784444
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.08160254037844385,
"y": -0.05866025403784443
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.09459292143521042,
"y": -0.06616025403784444
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.09459292143521042,
"y": -0.06616025403784444
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.09125052303291258,
"y": -0.07781549419036826
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.09459292143521042,
"y": -0.06616025403784444
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.09125052303291258,
"y": -0.07781549419036826
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.09125052303291258,
"y": -0.07781549419036826
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.12957057706115918,
"y": -0.11613554821861488
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.09125052303291258,
"y": -0.07781549419036826
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.12957057706115918,
"y": -0.11613554821861488
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.12957057706115918,
"y": -0.11613554821861488
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.16536162357033266,
"y": -0.05414363720866581
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.1506884202584923,
"y": -0.08700000000000009
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 7.005207871384369e-16,
"y": 4.184373959199281e-16
},
"radius": 0.17400000000000082,
"is_trigo": False,
"angle1": -2.825168403061107,
"angle2": -2.4108193529218807,
"angle": 0.4143490501392262
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.16536162357033266,
"y": -0.05414363720866581
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.11301545629335556,
"y": -0.040117523959935314
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.16536162357033266,
"y": -0.05414363720866581
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.11301545629335556,
"y": -0.040117523959935314
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.11301545629335556,
"y": -0.040117523959935314
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.10459292143521043,
"y": -0.048839745962155665
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.11301545629335556,
"y": -0.040117523959935314
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.10459292143521043,
"y": -0.048839745962155665
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.10459292143521043,
"y": -0.048839745962155665
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.09160254037844386,
"y": -0.04133974596215566
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.10459292143521043,
"y": -0.048839745962155665
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.09160254037844386,
"y": -0.04133974596215566
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.09160254037844386,
"y": -0.04133974596215566
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.10000000000000002,
"y": -0.010000000000000064
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.09400439217788162,
"y": -0.035541162773142815
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -6.442408349218697e-16,
"y": -1.172716245458817e-16
},
"radius": 0.10049875621120828,
"is_trigo": False,
"angle1": -3.0419240010986313,
"angle2": -2.7176625304826545,
"angle": 0.32426147061597677
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.1,
"y": -0.01000000000000002
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.115,
"y": -0.010000000000000021
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.1,
"y": -0.01000000000000002
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.115,
"y": -0.010000000000000021
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.115,
"y": -0.010000000000000021
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.11793301815030346,
"y": -0.02176493326044298
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.115,
"y": -0.010000000000000021
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.11793301815030346,
"y": -0.02176493326044298
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.11793301815030346,
"y": -0.02176493326044298
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.17027918542728054,
"y": -0.03579104650917346
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.11793301815030346,
"y": -0.02176493326044298
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.17027918542728054,
"y": -0.03579104650917346
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.17027918542728054,
"y": -0.03579104650917346
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.17027918542728054,
"y": 0.03579104650917339
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.17400000000000002,
"y": -3.196328145774592e-17
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 2.1854129343774616e-16,
"y": -2.2719414732835864e-17
},
"radius": 0.17400000000000024,
"is_trigo": False,
"angle1": 2.9344181285201794,
"angle2": -2.9344181285201794,
"angle": 0.4143490501392275
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.17027918542728054,
"y": 0.03579104650917339
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.11793301815030346,
"y": 0.021764933260442938
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.17027918542728054,
"y": 0.03579104650917339
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.11793301815030346,
"y": 0.021764933260442938
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.11793301815030346,
"y": 0.021764933260442938
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.115,
"y": 0.00999999999999998
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.11793301815030346,
"y": 0.021764933260442938
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.115,
"y": 0.00999999999999998
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.115,
"y": 0.00999999999999998
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.1,
"y": 0.009999999999999981
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.115,
"y": 0.00999999999999998
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.1,
"y": 0.009999999999999981
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.1,
"y": 0.009999999999999981
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.0916025403784439,
"y": 0.041339745962155595
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.09918077307993205,
"y": 0.01622264624736137
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": 5.149715067397098e-16,
"y": -8.630974361485694e-17
},
"radius": 0.10049875621120943,
"is_trigo": False,
"angle1": 2.717662530482658,
"angle2": 3.0419240010986313,
"angle": 0.3242614706159732
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.0916025403784439,
"y": 0.04133974596215555
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.10459292143521048,
"y": 0.04883974596215554
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.0916025403784439,
"y": 0.04133974596215555
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.10459292143521048,
"y": 0.04883974596215554
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.10459292143521048,
"y": 0.04883974596215554
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.1130154562933556,
"y": 0.040117523959935175
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.10459292143521048,
"y": 0.04883974596215554
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.1130154562933556,
"y": 0.040117523959935175
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.1130154562933556,
"y": 0.040117523959935175
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.16536162357033268,
"y": 0.05414363720866561
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.1130154562933556,
"y": 0.040117523959935175
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.16536162357033268,
"y": 0.05414363720866561
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.16536162357033268,
"y": 0.05414363720866561
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.12957057706115932,
"y": 0.11613554821861471
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.1506884202584924,
"y": 0.0869999999999999
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -1.8976277816168484e-16,
"y": 1.3754265283906457e-16
},
"radius": 0.1739999999999998,
"is_trigo": False,
"angle1": 2.410819352921881,
"angle2": 2.8251684030611095,
"angle": 0.4143490501392284
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.12957057706115932,
"y": 0.11613554821861471
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.09125052303291267,
"y": 0.07781549419036814
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.12957057706115932,
"y": 0.11613554821861471
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.09125052303291267,
"y": 0.07781549419036814
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.09125052303291267,
"y": 0.07781549419036814
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.0945929214352105,
"y": 0.06616025403784431
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.09125052303291267,
"y": 0.07781549419036814
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.0945929214352105,
"y": 0.06616025403784431
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.0945929214352105,
"y": 0.06616025403784431
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.08160254037844392,
"y": 0.05866025403784432
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.0945929214352105,
"y": 0.06616025403784431
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.08160254037844392,
"y": 0.05866025403784432
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.08160254037844392,
"y": 0.05866025403784432
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.05866025403784446,
"y": 0.08160254037844383
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.0777817459305203,
"y": 0.06363961030678922
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -1.6279046601097631e-15,
"y": 1.6846140302077266e-15
},
"radius": 0.1004987562112066,
"is_trigo": False,
"angle1": 2.194063754884355,
"angle2": 2.5183252255003374,
"angle": 0.32426147061598254
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.05866025403784443,
"y": 0.08160254037844383
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.06616025403784444,
"y": 0.09459292143521042
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.05866025403784443,
"y": 0.08160254037844383
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.06616025403784444,
"y": 0.09459292143521042
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.06616025403784444,
"y": 0.09459292143521042
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.07781549419036826,
"y": 0.09125052303291256
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.06616025403784444,
"y": 0.09459292143521042
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.07781549419036826,
"y": 0.09125052303291256
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.07781549419036826,
"y": 0.09125052303291256
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.11613554821861487,
"y": 0.12957057706115915
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.07781549419036826,
"y": 0.09125052303291256
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.11613554821861487,
"y": 0.12957057706115915
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.11613554821861487,
"y": 0.12957057706115915
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.054143637208665815,
"y": 0.16536162357033263
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.08700000000000009,
"y": 0.1506884202584923
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -1.9113612032239482e-16,
"y": 4.279936410867517e-16
},
"radius": 0.17399999999999957,
"is_trigo": False,
"angle1": 1.8872205773235815,
"angle2": 2.301569627462811,
"angle": 0.4143490501392293
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.054143637208665815,
"y": 0.16536162357033263
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.04011752395993532,
"y": 0.11301545629335555
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.054143637208665815,
"y": 0.16536162357033263
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.04011752395993532,
"y": 0.11301545629335555
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.04011752395993532,
"y": 0.11301545629335555
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.048839745962155665,
"y": 0.10459292143521043
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.04011752395993532,
"y": 0.11301545629335555
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.048839745962155665,
"y": 0.10459292143521043
}
},
{
"name": "",
"object_class": "volmdlr.edges.LineSegment2D",
"package_version": "0.1.13.dev167",
"points": [
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.048839745962155665,
"y": 0.10459292143521043
},
{
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.04133974596215566,
"y": 0.09160254037844384
}
],
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.048839745962155665,
"y": 0.10459292143521043
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.04133974596215566,
"y": 0.09160254037844384
}
},
{
"name": "",
"object_class": "volmdlr.edges.Arc2D",
"package_version": "0.1.13.dev167",
"start": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.04133974596215566,
"y": 0.09160254037844384
},
"end": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.010000000000000064,
"y": 0.1
},
"interior": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -0.035541162773142815,
"y": 0.09400439217788162
},
"center": {
"name": "",
"object_class": "volmdlr.core_compiled.Point2D",
"package_version": "0.1.13.dev167",
"x": -3.9129794287779357e-16,
"y": 1.6669209775476033e-15
},
"radius": 0.10049875621120721,
"is_trigo": False,
"angle1": 1.670464979286057,
"angle2": 1.9947264499020372,
"angle": 0.3242614706159801
}
]
})

rotor_internal=rotor_internal.primitives[0]
stator_external=stator_external.primitives[0]


rotor_mesher=vmmesh.Mesher([rotor_internal],[rotor_magnet,rotor_external],60)

surface=faces.Surface2D(stator_external,[stator_internal])
pattern=surface.get_pattern_single_inner()
pattern_mesher=vmmesh.Mesher([],[pattern],60)
all_patterns=surface.contour_from_pattern()
stator_mesher=vmmesh.Mesher([],all_patterns,60)
