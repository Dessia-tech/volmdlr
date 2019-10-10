#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 10:57:51 2019

@author: ringhausen
"""

import volmdlr as  vm

moteur = vm.Step('/home/ringhausen/Documents/git/ClientsProjects/Renault/CMO/data/step/MOTEUR HRevoUNIFY.stp')
boite = vm.Step('/home/ringhausen/Documents/git/ClientsProjects/Renault/CMO/data/step/BOITE E-TECHg2.stp')
cmo = vm.Step('/home/ringhausen/Documents/git/ClientsProjects/Renault/CMO/data/step/CMO2.stp')

volumemodel = cmo.to_volume_model('cmo')
boite.to_volume_model('boite', volumemodel)
moteur.to_volume_model('moteur', volumemodel)

shell0 = volumemodel.shells[0]
shell1 = volumemodel.shells[1]
shell2 = volumemodel.shells[2]

#pc = shell0.shell_intersection(shell1)
#print('intersection à', pc*100, '%')

shell0.Translation((0.05,0.1,0.3))
shell1.Translation((-0.1,-0.1,-0.1))
mesure = shell0.distance_to_shell(shell1, volumemodel)

volumemodel.BabylonShow()