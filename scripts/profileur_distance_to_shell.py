#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 10:57:51 2019

@author: ringhausen
"""

import volmdlr as  vm

moteur = vm.Step('/home/ringhausen/Documents/git/ClientsProjects/Renault/CMO/data/step/MOTEUR HRevoUNIFY.stp')
boite = vm.Step('/home/ringhausen/Documents/git/ClientsProjects/Renault/CMO/data/step/BOITE E-TECHg2.stp')
#cmo = vm.Step('/home/ringhausen/Documents/git/ClientsProjects/Renault/CMO/data/step/CMO2.stp')

volumemodel = moteur.to_volume_model('moteur')
boite.to_volume_model('boite', volumemodel)
#cmo.to_volume_model('cmo', volumemodel)

shell0 = volumemodel.shells[0] # LE MOTEUR
shell1 = volumemodel.shells[1] # LA BOITE
#shell2 = volumemodel.shells[2] # LE CMO
#shell3 = volumemodel.shells[3] # UN PETIT BOUT SE TROUVANT A L'INTERIEUR DU CMO
#shell4 = volumemodel.shells[4] # UN PETIT BOUT SE TROUVANT A L'INTERIEUR DU CMO

#del volumemodel.shells[4]
#del volumemodel.shells[3]

#pc = shell0.shell_intersection(shell1)
#print('intersection à', pc*100, '%')

shell0.Translation((0.05,0.1,0.6), False)
#shell1.Translation((-0.1,-0.1,-0.1))
#mesure = shell0.distance_to_shell(shell1, volumemodel)
#print(mesure.distance)
#volumemodel.BabylonShow()