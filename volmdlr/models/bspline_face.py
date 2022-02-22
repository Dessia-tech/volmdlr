#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: s.bendjebla
"""

# %% Librairies

import geomdl
import volmdlr as vm
import volmdlr.faces as vmf

# %% Control points

control_points = [vm.Point3D(2.241, 0, 1.094),
                  vm.Point3D(2.241, 0.113, 1.094),
                  vm.Point3D(2.241, 0.223, 1.082),
                  vm.Point3D(2.241, 0.330, 1.059),
                  vm.Point3D(2.241, 0.536, 0.990),
                  vm.Point3D(2.241, 0.723, 0.879),
                  vm.Point3D(2.241, 0.811, 0.814),
                  vm.Point3D(2.241, 0.976, 0.665),
                  vm.Point3D(2.241, 1.112, 0.482),
                  vm.Point3D(2.241, 1.172, 0.383),
                  vm.Point3D(2.241, 1.274, 0.171),
                  vm.Point3D(2.241, 1.337, -0.064),
                  vm.Point3D(2.241, 1.359, -0.188),
                  vm.Point3D(2.241, 1.370, -0.317),
                  vm.Point3D(2.241, 1.370, -0.452), 
                  
                  vm.Point3D(2.593, 0, 1.323),
                  vm.Point3D(2.593, 0.120, 1.323),
                  vm.Point3D(2.592, 0.240, 1.312),
                  vm.Point3D(2.592, 0.358, 1.289),
                  vm.Point3D(2.591, 0.587, 1.218),
                  vm.Point3D(2.592, 0.798, 1.098),
                  vm.Point3D(2.594, 0.897, 1.0264),
                  vm.Point3D(2.598, 1.081, 0.860),
                  vm.Point3D(2.605, 1.231, 0.656),
                  vm.Point3D(2.609, 1.295, 0.545),
                  vm.Point3D(2.617, 1.402, 0.309),
                  vm.Point3D(2.624, 1.466, 0.0483),
                  vm.Point3D(2.628, 1.486, -0.086),
                  vm.Point3D(2.630, 1.496, -0.226),
                  vm.Point3D(2.632, 1.497, -0.369),
                            
                  vm.Point3D(2.969, 0, 1.506),
                  vm.Point3D(2.969, 0.127, 1.506),
                  vm.Point3D(2.968, 0.254, 1.495),
                  vm.Point3D(2.967, 0.380, 1.472),
                  vm.Point3D(2.967, 0.628, 1.398),
                  vm.Point3D(2.969, 0.858, 1.271),
                  vm.Point3D(2.971, 0.966, 1.194),
                  vm.Point3D(2.978, 1.166, 1.016),
                  vm.Point3D(2.988, 1.327, 0.796),
                  vm.Point3D(2.994, 1.397, 0.676),
                  vm.Point3D(3.007, 1.510, 0.422),
                  vm.Point3D(3.019, 1.575, 0.144),
                  vm.Point3D(3.024, 1.596, 0.001),
                  vm.Point3D(3.028, 1.606, -0.144),
                  vm.Point3D(3.030, 1.606, -0.293),  
                  
                  vm.Point3D(3.369, 0, 1.645),
                  vm.Point3D(3.369, 0.132, 1.645),
                  vm.Point3D(3.368, 0.265, 1.634),
                  vm.Point3D(3.367, 0.398, 1.610),
                  vm.Point3D(3.366, 0.660, 1.533),
                  vm.Point3D(3.369, 0.905, 1.401),
                  vm.Point3D(3.371, 1.020, 1.321),
                  vm.Point3D(3.378, 1.233, 1.135),
                  vm.Point3D(3.390, 1.404, 0.904),
                  vm.Point3D(3.396, 1.477, 0.779),
                  vm.Point3D(3.410, 1.597, 0.513),
                  vm.Point3D(3.423, 1.666, 0.224),
                  vm.Point3D(3.429, 1.687, 0.076),
                  vm.Point3D(3.433, 1.697, -0.073),
                  vm.Point3D(3.436, 1.698, -0.225),
                  
                  vm.Point3D(4.202, 0, 1.849),
                  vm.Point3D(4.202, 0.141, 1.849),
                  vm.Point3D(4.201, 0.284, 1.837),
                  vm.Point3D(4.200, 0.428, 1.811),
                  vm.Point3D(4.199, 0.711, 1.730),
                  vm.Point3D(4.201, 0.977, 1.590),
                  vm.Point3D(4.203, 1.102, 1.505),
                  vm.Point3D(4.209, 1.335, 1.309),
                  vm.Point3D(4.218, 1.522, 1.066),
                  vm.Point3D(4.224, 1.603, 0.934),
                  vm.Point3D(4.235, 1.735, 0.655),
                  vm.Point3D(4.246, 1.811, 0.354),
                  vm.Point3D(4.251, 1.835, 0.200),
                  vm.Point3D(4.255, 1.847, 0.044),
                  vm.Point3D(4.257, 1.847, -0.109),
                  
                  vm.Point3D(5.068, 0, 1.938),
                  vm.Point3D(5.068, 0.147, 1.938),
                  vm.Point3D(5.068, 0.296, 1.925),
                  vm.Point3D(5.067, 0.446, 1.897),
                  vm.Point3D(5.067, 0.739, 1.812),
                  vm.Point3D(5.068, 1.013, 1.668),
                  vm.Point3D(5.068, 1.143, 1.582),
                  vm.Point3D(5.071, 1.384, 1.384),
                  vm.Point3D(5.075, 1.580, 1.139),
                  vm.Point3D(5.077, 1.665, 1.007),
                  vm.Point3D(5.083, 1.807, 0.728),
                  vm.Point3D(5.087, 1.892, 0.426),
                  vm.Point3D(5.090, 1.919, 0.272),
                  vm.Point3D(5.091, 1.933, 0.117),
                  vm.Point3D(5.092, 1.933, -0.036),
                  
                  vm.Point3D(5.502, 0, 1.962),
                  vm.Point3D(5.502, 0.150, 1.962),
                  vm.Point3D(5.502, 0.301, 1.948),
                  vm.Point3D(5.502, 0.452, 1.920),
                  vm.Point3D(5.502, 0.749, 1.832),
                  vm.Point3D(5.502, 1.025, 1.687),
                  vm.Point3D(5.502, 1.155, 1.601),
                  vm.Point3D(5.503, 1.398, 1.402),
                  vm.Point3D(5.505, 1.597, 1.159),
                  vm.Point3D(5.506, 1.683, 1.028),
                  vm.Point3D(5.508, 1.829, 0.750),
                  vm.Point3D(5.510, 1.917, 0.449),
                  vm.Point3D(5.511, 1.946, 0.296),
                  vm.Point3D(5.511, 1.961, 0.141),
                  vm.Point3D(5.512, 1.961, -0.012),
                  
                  vm.Point3D(5.929, 0, 1.974),
                  vm.Point3D(5.929, 0.152, 1.974),
                  vm.Point3D(5.929, 0.306, 1.960),
                  vm.Point3D(5.929, 0.458, 1.930),
                  vm.Point3D(5.930, 0.756, 1.841),
                  vm.Point3D(5.930, 1.032, 1.694),
                  vm.Point3D(5.930, 1.162, 1.607),
                  vm.Point3D(5.930, 1.405, 1.409),
                  vm.Point3D(5.930, 1.604, 1.167),
                  vm.Point3D(5.930, 1.691, 1.036),
                  vm.Point3D(5.931, 1.839, 0.760),
                  vm.Point3D(5.931, 1.929, 0.461),
                  vm.Point3D(5.931, 1.960, 0.308),
                  vm.Point3D(5.931, 1.975, 0.153),
                  vm.Point3D(5.9315, 1.975, 0),
                  
                  vm.Point3D(6.35, 0, 1.974),
                  vm.Point3D(6.35, 0.155, 1.974),
                  vm.Point3D(6.35, 0.310, 1.959),
                  vm.Point3D(6.35, 0.463, 1.929),
                  vm.Point3D(6.35, 0.761, 1.838),
                  vm.Point3D(6.35, 1.036, 1.691),
                  vm.Point3D(6.35, 1.166, 1.605),
                  vm.Point3D(6.35, 1.407, 1.407),
                  vm.Point3D(6.35, 1.605, 1.166),
                  vm.Point3D(6.35, 1.691, 1.036),
                  vm.Point3D(6.35, 1.838, 0.761),
                  vm.Point3D(6.35, 1.9293, 0.463),
                  vm.Point3D(6.35, 1.959, 0.310),
                  vm.Point3D(6.35, 1.975, 0.155),
                  vm.Point3D(6.35, 1.975, 0)] 

# %% Bspline-surface parameters 

degree_u, degree_v, nb_u, nb_v = 5, 5, 9, 15
knots_vector_u = geomdl.knotvector.generate(degree_u, nb_u)
knots_vector_v = geomdl.knotvector.generate(degree_v, nb_v)

(u_knots, u_multiplicities) = vmf.knots_vector_inv(knots_vector_u)
(v_knots, v_multiplicities) = vmf.knots_vector_inv(knots_vector_v)

# %% Bspline-surface definition

bspline_surface = vm.faces.BSplineSurface3D(degree_u = degree_u, 
                                            degree_v = degree_v, 
                                            control_points = control_points,
                                            nb_u = nb_u, 
                                            nb_v = nb_v,
                                            u_multiplicities = u_multiplicities, 
                                            v_multiplicities = v_multiplicities, 
                                            u_knots = u_knots, 
                                            v_knots = v_knots)

# %% Display 

# bspline_surface.plot()

bspline_face = bspline_surface.rectangular_cut(0, 1, 0, 1)
# bspline_face.babylonjs()

