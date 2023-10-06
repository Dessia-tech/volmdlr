#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 2022.

@author: s.bendjebla
"""

# %% Libraries
import volmdlr
from volmdlr import edges, wires

# %% Contour2d_1

primitives = [
    edges.LineSegment2D(volmdlr.Point2D(0.001, 0.014),
                        volmdlr.Point2D(0.001, 0.0125)),

    edges.Arc2D.from_3_points(volmdlr.Point2D(0.001, 0.0125),
                              volmdlr.Point2D(0.009862829911410362, 0.007744326060968065),
                              volmdlr.Point2D(0.012539936203984454, 0.0)),

    edges.Arc2D.from_3_points(volmdlr.Point2D(0.012539936203984454, 0.0),
                              volmdlr.Point2D(0.0, -0.012539936203984454),
                              volmdlr.Point2D(-0.012539936203984454, 0.0)),

    edges.Arc2D.from_3_points(volmdlr.Point2D(-0.012539936203984454, 0.0),
                              volmdlr.Point2D(-0.00921384654213387, 0.008506176103162205),
                              volmdlr.Point2D(-0.001, 0.0125)),

    edges.LineSegment2D(volmdlr.Point2D(-0.001, 0.0125),
                        volmdlr.Point2D(-0.001, 0.014)),

    edges.LineSegment2D(volmdlr.Point2D(-0.001, 0.014),
                        volmdlr.Point2D(0.001, 0.014))
]

contour2d_1 = wires.Contour2D(primitives)

# %% Contour2d_2

points2d = [volmdlr.Point2D(0.7557261768236382, 0.047840610095316816),
            volmdlr.Point2D(0.7557153778090784, 0.053529840943971466),
            volmdlr.Point2D(0.7548352148318592, 0.059379242604894106),
            volmdlr.Point2D(0.7509972673886837, 0.06812701209830616),
            volmdlr.Point2D(0.7463815515547884, 0.06876367708780329),
            volmdlr.Point2D(0.7426873841819802, 0.06850598841747624),
            volmdlr.Point2D(0.7386859822441963, 0.0690723949917957),
            volmdlr.Point2D(0.7344087842895368, 0.0660771254673235),
            volmdlr.Point2D(0.7316087522378829, 0.0572133737658937),
            volmdlr.Point2D(0.7306755405178993, 0.04995246337945117),
            volmdlr.Point2D(0.7310451707836978, 0.044346281865748856)]

curve2d_1 = edges.BSplineCurve2D(degree=3,
                                 control_points=points2d,
                                 knot_multiplicities=[4, 1, 1, 1, 1, 1, 1, 1, 4],
                                 knots=[0.0, 0.26206414743620277, 0.3531225286567413,
                                        0.42241139891509305, 0.4822254408401885,
                                        0.5440692281614715, 0.6202023477514225,
                                        0.7222748100879735, 1.0])

points2d = [volmdlr.Point2D(0.7310451707836978, 0.044346281865748856),
            volmdlr.Point2D(0.7308827017192874, 0.03840965331596735),
            volmdlr.Point2D(0.7318250470334302, 0.03228668394110015),
            volmdlr.Point2D(0.7354633489443133, 0.023006441447214778),
            volmdlr.Point2D(0.7402551819180077, 0.022102466430262133),
            volmdlr.Point2D(0.7440795569187786, 0.02242177670352842),
            volmdlr.Point2D(0.7482500898246996, 0.02178574616551938),
            volmdlr.Point2D(0.7525849938211369, 0.025193737559804485),
            volmdlr.Point2D(0.7553036624157348, 0.03437896214237193),
            volmdlr.Point2D(0.7561315600530795, 0.042073178367651974),
            volmdlr.Point2D(0.7557261768236382, 0.047840610095316816)]

curve2d_2 = edges.BSplineCurve2D(degree=3,
                                 control_points=points2d,
                                 knot_multiplicities=[4, 1, 1, 1, 1, 1, 1, 1, 4],
                                 knots=[0.0, 0.26231860269644847, 0.35446984209379995,
                                        0.4237193471508679, 0.48275286541966095,
                                        0.5439383602451768, 0.62045686469776,
                                        0.7233434357617774, 1.0])

contour2d_2 = wires.Contour2D([curve2d_1, curve2d_2])

points = [
    volmdlr.Point2D(0.20308817713481986, 0.04966773764193705),
    volmdlr.Point2D(0.7969119765656952, 0.04966797879396336),
    volmdlr.Point2D(0.8442697800221348, 0.04966805448106126),
    volmdlr.Point2D(0.9031379981759905, 0.04545485716230839),
    volmdlr.Point2D(0.9479619187253645, 0.10517762923653938),
    volmdlr.Point2D(0.9545454713622077, 0.19659816573695035),
    volmdlr.Point2D(0.9193352753725873, 0.27036785408838365),
    volmdlr.Point2D(0.9193352733808331, 0.2891129034169921),
    volmdlr.Point2D(0.9193352491395294, 0.7108874440320542),
    volmdlr.Point2D(0.9193352445733572, 0.7296301344365137),
    volmdlr.Point2D(0.9545453968718225, 0.8034021061594617),
    volmdlr.Point2D(0.9479618025250582, 0.8948226395224322),
    volmdlr.Point2D(0.9031378875027738, 0.9545453756361115),
    volmdlr.Point2D(0.8442696587355882, 0.9503323074053969),
    volmdlr.Point2D(0.7969118264212707, 0.950332267803071),
    volmdlr.Point2D(0.203088054423346, 0.9503320443667395),
    volmdlr.Point2D(0.15573023404235076, 0.9503319931770652),
    volmdlr.Point2D(0.09686198756850553, 0.9545451189134176),
    volmdlr.Point2D(0.05203811911069425, 0.8948223672618099),
    volmdlr.Point2D(0.04545450905371345, 0.8034018054641981),
    volmdlr.Point2D(0.08066468979315351, 0.7296322065226041),
    volmdlr.Point2D(0.08066468965825305, 0.7108871704475791),
    volmdlr.Point2D(0.08066474955300275, 0.289112539074918),
    volmdlr.Point2D(0.08066475417207301, 0.27036750263798714),
    volmdlr.Point2D(0.04545459207437027, 0.19659788603626316),
    volmdlr.Point2D(0.05203818215846509, 0.10517733497317228),
    volmdlr.Point2D(0.09686210943003962, 0.04545460320088811),
    volmdlr.Point2D(0.15573034610061637, 0.04966777166777087)]

contour1_cut_by_wire = wires.Contour2D.from_points(points)

# %% Contour 2

points = [
    volmdlr.Point2D(0.2030881575366132, 0.04966771677601732),
    volmdlr.Point2D(0.20308809125575447, 0.2891126765655333),
    volmdlr.Point2D(0.08066474910267005, 0.2891125988482983),
    volmdlr.Point2D(0.05674291581332371, 0.28911260559053886),
    volmdlr.Point2D(0.04267127288273421, 0.311080039363985),
    volmdlr.Point2D(0.04267123136672513, 0.6889196955746452),
    volmdlr.Point2D(0.056742861492621116, 0.7108870816261534),
    volmdlr.Point2D(0.08066469250960157, 0.7108871106086502),
    volmdlr.Point2D(0.2030880435632797, 0.7108871156873701),
    volmdlr.Point2D(0.20308804388447252, 0.9503320678576602),
    volmdlr.Point2D(0.20308798646643267, 0.9876766506444887),
    volmdlr.Point2D(0.21715969201684732, 1.0),
    volmdlr.Point2D(0.7828401610275215, 1.0),
    volmdlr.Point2D(0.7969118753894622, 0.9876767973117213),
    volmdlr.Point2D(0.7969118376396574, 0.9503322697795179),
    volmdlr.Point2D(0.796911890394745, 0.7108872963268008),
    volmdlr.Point2D(0.9193352406490979, 0.7108873750750229),
    volmdlr.Point2D(0.943257068363107, 0.7108873714262228),
    volmdlr.Point2D(0.957328707411281, 0.6889199122213159),
    volmdlr.Point2D(0.9573287795940568, 0.3110803289147694),
    volmdlr.Point2D(0.9432571269950443, 0.28911295065634096),
    volmdlr.Point2D(0.9193352813560651, 0.2891129208514467),
    volmdlr.Point2D(0.7969119478569466, 0.289112868790721),
    volmdlr.Point2D(0.7969119771240915, 0.049667935920999225),
    volmdlr.Point2D(0.7969119981329132, 0.012323366834239619),
    volmdlr.Point2D(0.7828402953854121, 0.0),
    volmdlr.Point2D(0.21715982313378532, 0.0),
    volmdlr.Point2D(0.2030881157414971, 0.012323183603131671)]

contour2_cut_by_wire = wires.Contour2D.from_points(points)

line_segment1 = edges.LineSegment2D(volmdlr.Point2D(1, -1), volmdlr.Point2D(1.5, 1))
arc = edges.Arc2D.from_3_points(volmdlr.Point2D(1.5, 1), volmdlr.Point2D(1.3, 1.5), volmdlr.Point2D(0.5, 1.5))
line_segment2 = edges.LineSegment2D(volmdlr.Point2D(0.5, 1.5), volmdlr.Point2D(-2, 1))
line_segment3 = edges.LineSegment2D(volmdlr.Point2D(-2, 1), volmdlr.Point2D(-2, 0.7))
line_segment4 = edges.LineSegment2D(volmdlr.Point2D(-2, 0.7), volmdlr.Point2D(-1, 1))
points2d = [volmdlr.Point2D(-1, 1), volmdlr.Point2D(2, 2), volmdlr.Point2D(-2, -2), volmdlr.Point2D(1, -1)]
bspline = edges.BSplineCurve2D(3, points2d, knot_multiplicities=[4, 4], knots=[0.0, 1.0])
contour2_unittest = wires.Contour2D([bspline, line_segment1, arc, line_segment2, line_segment3, line_segment4])
unordered_contour2_unittest = wires.Contour2D([line_segment2, bspline.reverse(), arc.reverse(),
                                                          line_segment1, line_segment3, line_segment4])

invalid_unordered_contour2_unittest = wires.Contour2D([line_segment2, bspline.reverse(), arc.reverse(),
                                                       line_segment1, line_segment3, line_segment4,
                                                       edges.LineSegment2D(volmdlr.Point2D(1, -1),
                                                                           volmdlr.Point2D(1.5, -1))])

unordered_wire_unittest = wires.Wire2D([line_segment2, bspline.reverse(), arc.reverse(),
                                        line_segment1, line_segment3])

# Contour3D with all edges.
control_points = [volmdlr.Point3D(0, 3, 0),
                  volmdlr.Point3D(3, 2, 1),
                  volmdlr.Point3D(5, -1, 4),
                  volmdlr.Point3D(5, -4, 0),
                  volmdlr.Point3D(-1, -2, -3),
                  volmdlr.Point3D(-3, 4, 1)]
knots = [0.0, 1.0]
knot_multiplicities = [6, 6]
bspline_curve3d = edges.BSplineCurve3D(degree=5, control_points=control_points,
                                       knot_multiplicities=knot_multiplicities,
                                       knots=knots,
                                       weights=None,
                                       name='B Spline Curve 3D 1')
lineseg1 = edges.LineSegment3D(volmdlr.Point3D(3, 3, 2), bspline_curve3d.start)
lineseg2 = edges.LineSegment3D(bspline_curve3d.end, volmdlr.Point3D(-3, -3, 0))
arc = edges.Arc3D.from_3_points(volmdlr.Point3D(-3, -3, 0),
                                volmdlr.Point3D(6.324555320336761, -5.692099788303083, -0.8973665961010275),
                                volmdlr.Point3D(3, 3, 2))
contour3d = wires.Contour3D([lineseg1, bspline_curve3d, lineseg2, arc])
