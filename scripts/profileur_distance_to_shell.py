#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 10:57:51 2019

@author: ringhausen
"""

import volmdlr as  vm
import volmdlr.primitives3D as p3d

moteur = vm.Step('/home/ringhausen/Documents/git/ClientsProjects/Renault/CMO/data/step/MOTEUR HRevoUNIFY.stp')
boite = vm.Step('/home/ringhausen/Documents/git/ClientsProjects/Renault/CMO/data/step/BOITE E-TECHg2.stp')
#cmo = vm.Step('/home/ringhausen/Documents/git/ClientsProjects/Renault/CMO/data/step/CMO2.stp')

volumemodel = moteur.to_volume_model('moteur')
boite.to_volume_model('boite', volumemodel)
#cmo.to_volume_model('cmo', volumemodel)

#origin = vm.Point3D((0,0,0))
#u = vm.Vector3D((2,1,0))
#v = vm.Vector3D((0,1,0))
#w = vm.Vector3D((0.5,0,1))
#frame = vm.Frame3D(origin, u, v, w)
#block = p3d.Block(frame)
#volumemodel.shells.append(block)


shell0 = volumemodel.shells[0] # LE MOTEUR
shell1 = volumemodel.shells[1] # LA BOITE
#shell2 = volumemodel.shells[2] # LE BLOCK
#shell2 = volumemodel.shells[2] # LE CMO
#shell3 = volumemodel.shells[3] # UN PETIT BOUT SE TROUVANT A L'INTERIEUR DU CMO
#shell4 = volumemodel.shells[4] # UN PETIT BOUT SE TROUVANT A L'INTERIEUR DU CMO

#del volumemodel.shells[4]
#del volumemodel.shells[3]


#%%

origin = vm.Point3D((1,0.1,0.6))
u = vm.y3D
v = vm.z3D
w = vm.x3D
frame = vm.Frame3D(origin, u, v, w)



shell0.frame_mapping(frame, 'new', False)
shell1.Translation((0.5,0,0), False)

volumemodel.shells[0] = shell0.frame_mapping(frame, 'new', True)
volumemodel.shells[1] = shell1.Translation((0.5,0,0), True)

#volumemodel.shells.append(volumemodel.shells[0].bounding_box)
#volumemodel.shells.append(volumemodel.shells[1].bounding_box)

#for face in volumemodel.shells[0].faces+volumemodel.shells[1].faces:
#    volumemodel.shells.append(face.bounding_box)



mesure = volumemodel.shells[0].distance_to_shell(volumemodel.shells[1], volumemodel)
if mesure is not None:
    volumemodel.shells.append(mesure)
    print(mesure.distance)

#pc = shell0.shell_intersection(shell1)
#print('intersection à', pc*100, '%')

volumemodel.BabylonShow()

#%%