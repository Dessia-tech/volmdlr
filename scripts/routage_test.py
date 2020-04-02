#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 16:55:25 2019

@author: ringhausen
"""

import volmdlr as  vm
import volmdlr.primitives3D as p3d


### BLOCK DE DEPART ###
origin = vm.Point3D((-2,1,0))
u = vm.Vector3D((0.2,0,0))
v = vm.Vector3D((0,0.2,0))
w = vm.Vector3D((0,0,0.2))
frame = vm.Frame3D(origin, u, v, w)
block_depart = p3d.Block(frame, color=(0,0.5,0))


### BLOCK D'ARIVEE ###
origin = vm.Point3D((2,-1,0))
u = vm.Vector3D((0.2,0,0))
v = vm.Vector3D((0,0.2,0))
w = vm.Vector3D((0,0,0.2))
frame = vm.Frame3D(origin, u, v, w)
block_arrivee = p3d.Block(frame, color=(0,0,0.5))


### BLOCKS OBSTACLES ###
origin = vm.Point3D((-0.5,-0.4,0))
u = vm.Vector3D((1,0,0))
v = vm.Vector3D((0,0.9,0))
w = vm.Vector3D((0,0,0.6))
frame = vm.Frame3D(origin, u, v, w)
block_obstacle1 = p3d.Block(frame, color=(0.5,0,0))

origin = vm.Point3D((0.7,-0.5,0))
u = vm.Vector3D((3,0,0))
v = vm.Vector3D((0,0.5,0))
w = vm.Vector3D((0,0,2))
frame = vm.Frame3D(origin, u, v, w)
block_obstacle2 = p3d.Block(frame, color=(0.5,0,0))


### POINT DE DEPART ET D'ARRIVEE DU ROUTAGE
x1 = block_depart.frame.origin[0] + block_depart.frame.u[0]/2
y1 = block_depart.frame.origin[1] - 0.8
z1 = block_depart.frame.origin[2]
point_depart = vm.Point3D((x1, y1, z1))

x2 = block_arrivee.frame.origin[0] - block_arrivee.frame.u[0]/2
y2 = block_arrivee.frame.origin[1] 
z2 = block_arrivee.frame.origin[2]+0.2
point_arrivee = vm.Point3D((x2, y2, z2))


# ATTENTION : nouvelle façon de definir un VolumeModel
# L'attribut primitives a regoupe les Shell3D, les Mesure, etc.
# On peut tout de même appeler l'attribut VolumeModel.shells qui va chercher 
# dans VolumeModel.primitives les objets Shell3D
primitives = [block_depart, block_arrivee, block_obstacle1, block_obstacle2, point_depart, point_arrivee]
volumemodel = vm.VolumeModel(primitives, 'test')


### CREATION DU ROUTAGE ###
routage = vm.Routing(point_depart, point_arrivee, volumemodel=volumemodel)
### ALGORITHME 1 ###
#good_mesures, bad_mesures = routage.straight_line()
### AJOUT DANS LES PRIMITIVES POUR LE BABYLONSHOW ###
#volumemodel.primitives.extend(good_mesures)
#volumemodel.primitives.extend(bad_mesures)
#good_distance = sum([mesure.distance for mesure in good_mesures])
#bad_distance = sum([mesure.distance for mesure in bad_mesures])
#print('distance en bleu', good_distance)
#print('distance en rouge', bad_distance)
### ALGORITHME 2 ###
mesures = routage.straight_line()
### AJOUT DANS LES PRIMITIVES POUR LE BABYLONSHOW ###
volumemodel.primitives.extend(mesures)

volumemodel.babylonjs_from_script()