#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 10:57:51 2019

@author: ringhausen
"""

import volmdlr as  vm

moteur = vm.Step('/home/ringhausen/Bureau/Renault/MOTEUR.txt')
boite = vm.Step('/home/ringhausen/Bureau/Renault/BOITE.txt')

volumemodel = moteur.to_volume_model('moteur')
boite.to_volume_model('boite', volumemodel)

shell0 = volumemodel.shells[0]
shell1 = volumemodel.shells[1]

pc = shell0.shell_intersection(shell1)
print('intersection à', pc*100, '%')

#shell0.Translation((600,100,0))
#shell0.distance_to_shell(shell1, volumemodel)

