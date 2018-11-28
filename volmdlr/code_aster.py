#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 14:08:23 2017

@author: steven
"""
import numpy as npy
import os


class CylinderNodeGroup:
    def __init__(self, name, mesh, point, normal, radius, inside=True):
        self.point = point
        self.normal = normal
        self.name = name
        self.mesh = mesh
        if inside:
            self.radius = 0
            self.precision = radius
        else:
            self.radius = radius
            self.precision = radius/100.
            
    def CodeAsterExport(self, a=None):
        char = '''DEFI_GROUP( reuse  = {}, MAILLAGE = {},
            CREA_GROUP_NO = _F(NOM ='{}', OPTION = 'ENV_CYLINDRE',
            POINT = {}, VECT_NORMALE = {}, RAYON = {}, PRECISION = {},),); \n \n''' \
            .format(self.mesh.name, self.mesh.name, self.name, self.point, 
                    self.normal, self.radius, self.precision)
        if a is not None:
            a.write(char)
        else:
            return char
    
class PlanarNodeGroup:
    def __init__(self, name, mesh, point, normal):
        self.point = point
        self.normal = normal
        self.name = name
        self.mesh = mesh
        self.precision = 0.0001
        
    def CodeAsterExport(self, a=None):
        char = '''DEFI_GROUP( reuse  = {}, MAILLAGE = {},
            CREA_GROUP_NO = _F(NOM = '{}', OPTION = 'PLAN',
            POINT = {}, VECT_NORMALE = {}, PRECISION = {},),); \n \n''' \
            .format(self.mesh.name, self.mesh.name, self.name, self.point, 
                    self.normal, self.precision)
        if a is not None:
            a.write(char)
        else:
            return char
        
class NodeGroupAssembly:
    def __init__(self, name, mesh, typ, list_node_group):
        self.name = name
        self.mesh = mesh
        self.typ = typ
        self.list_node_group = list_node_group
        
    def CodeAsterExport(self, a=None):
        if self.typ == 'diff':
            char = '''DEFI_GROUP( reuse  = {}, MAILLAGE = {},
                CREA_GROUP_NO = _F(NOM = '{}', DIFFE = {},),); \n \n''' \
                .format(self.mesh.name, self.mesh.name, self.name, 
                        [n.name for n in self.list_node_group])
        elif self.typ == 'union':
            char = '''DEFI_GROUP( reuse  = {}, MAILLAGE = {}, 
                CREA_GROUP_NO = _F(NOM = '{}', UNION = {},),); \n \n''' \
                .format(self.mesh.name, self.mesh.name, self.name, 
                        [n.name for n in self.list_node_group])
        elif self.typ == 'intersec':
            char = '''DEFI_GROUP( reuse  = {}, MAILLAGE = {}, 
                CREA_GROUP_NO = _F(NOM = '{}', INTERSEC = {},),); \n \n''' \
                .format(self.mesh.name, self.mesh.name, self.name, 
                        [n.name for n in self.list_node_group])
        if a is not None:
            a.write(char)
        else:
            return char
        
class LimitedPlaneNodeGroup:
    def __init__(self, name, mesh, point, normal, list_dir, list_length):
        self.point = point
        self.normal = normal
        self.name = name
        self.mesh = mesh
        self.precision = 0.0001
        self.list_node_group = []
        list_intersec = []
        self.list_node_group.append(PlanarNodeGroup(name = self.name + '_0', mesh = self.mesh,
                                                  point = self.point, normal = self.normal))
        list_intersec.append(self.list_node_group[-1])
        compt = 1
        for i, li in enumerate(list_dir):
            list_union = []
            normal1 = tuple(npy.cross(normal, list_dir[i]))
            radius1 = abs(list_length[i][1]/2.)
            point1 = tuple(point + npy.array(list_dir[i])*list_length[i][1]/2.)
            self.list_node_group.append(CylinderNodeGroup(name = self.name + '_' + str(compt),
                                                          mesh = self.mesh, point = point1, 
                                                          normal = normal1, radius = radius1))
            compt += 1
            list_union.append(self.list_node_group[-1])
            normal2 = tuple(npy.cross(normal, list_dir[i]))
            radius2 = abs(list_length[i][0]/2.)
            point2 = tuple(point + npy.array(list_dir[i])*list_length[i][0]/2.)
            self.list_node_group.append(CylinderNodeGroup(name = self.name + '_' + str(compt),
                                                          mesh = self.mesh, point = point2, 
                                                          normal = normal2, radius = radius2))
            compt += 1
            list_union.append(self.list_node_group[-1])
            self.list_node_group.append(NodeGroupAssembly(name = self.name + '_' + str(compt),
                                                          mesh = self.mesh, typ = 'union', 
                                                          list_node_group = list_union))
            compt += 1
            list_intersec.append(self.list_node_group[-1])
        self.list_node_group.append(NodeGroupAssembly(name = self.name,
                                                      mesh = self.mesh, typ = 'intersec', 
                                                      list_node_group = list_intersec))
        
    def CodeAsterExport(self, a=None):
        if a is None:
            export = []
            for nd_group in self.list_node_group:
                export.append(nd_group.CodeAsterExport())
            return export
        else:
            for nd_group in self.list_node_group:
                nd_group.CodeAsterExport(a)
                
class Mesh:
    def __init__(self, name, file_name_CAD, file_name_mesh=None, 
                 gmsh_path='/Users/Pierrem/DessIA/PowerPack/scripts/Gmsh.app/Contents/MacOS/./gmsh'):
        if file_name_mesh is None:
            self.file_name_mesh = name + '.msh'
        else:
            self.file_name_mesh = file_name_mesh
        self.name = name
        self.file_name_CAD = file_name_CAD
        self.gmsh_path = gmsh_path
        self.GenereMesh()
        
    def CodeAsterExport(self, a=None):
        char = '''PRE_GMSH(); \n
{} = LIRE_MAILLAGE(FORMAT = "ASTER",); \n \n'''.format(self.name)
        if a is None:
            return char
        else:
            a.write(char)
            
    def GenereMesh(self):
        a = open(self.name + '.geo', 'w')
        
        a.write('Merge "{}"; \n'.format(self.file_name_CAD))
        a.write('bb() = BoundingBox Volume {1}; \n')
        a.write('BoundingBox {bb(0), bb(3), bb(1), bb(4), bb(2), bb(5)}; \n')
        a.write('Mesh.CharacteristicLengthFactor = 0.1; \n')
        a.write('Mesh.ScalingFactor=0.001; \n')
        a.write('Characteristic Length {:} = 0.15; \n')
        a.write('Save "{}"; \n'.format(self.file_name_mesh))
        
        a.close()
        
        arg = '{} -3 -o {} -order 2'.format(self.name + '.geo', self.file_name_mesh)
#        output=subprocess.call([self.gmsh_path, arg])
        os.system('{} {}'.format(self.gmsh_path, arg))
        

                
class Material:
    def __init__(self, name, E, nu, rho, mesh):
        self.name = name
        self.mesh = mesh
        self.E = E
        self.nu = nu
        self.rho = rho
        
    def CodeAsterExport(self, a=None):
        char = '''{} = DEFI_MATERIAU(ELAS = _F(E = {},
                          NU = {},
                          RHO = {},),); \n
{} = AFFE_MATERIAU(MAILLAGE = {},
                    AFFE = _F(TOUT = 'OUI',
                            MATER = {},),); \n \n'''.format(
                self.name + '_1', self.E, self.nu, self.rho, self.name, 
                self.mesh.name, self.name + '_1')
        if a is None:
            return char
        else:
            a.write(char)
        
class Model:
    def __init__(self, name, mesh):
        self.name = name
        self.mesh = mesh
    def CodeAsterExport(self, a=None):
        char = '''{} = AFFE_MODELE(MAILLAGE = {},
                   AFFE = _F(TOUT = 'OUI',
                           PHENOMENE = 'MECANIQUE',
                           MODELISATION = '3D',),); \n \n'''.format(
                           self.name, self.mesh.name)
        if a is None:
            return char
        else:
            a.write(char)
    
                
class BoundaryCondition:
    def __init__(self, name, group, modele, displacement=None, load=None):
        self.name = name
        self.group = group
        self.displacement = displacement
        self.modele = modele
        
    def CodeAsterExport(self, a=None):
        char_ddl = ''
        for D, val in self.displacement.items():
            if val is not None:
                char_ddl += D + ' = ' + str(val) + ', '
        char = '''{} = AFFE_CHAR_MECA(MODELE = {},
                       DDL_IMPO = _F(GROUP_NO = '{}',
                                   {}),); \n \n'''.format(
                   self.name, self.modele.name, self.group.name, char_ddl)
        if a is None:
            return char
        else:
            a.write(char)
            
class MecaStat:
    def __init__(self, name, modele, material, boundary):
        self.name = name
        self.modele = modele
        self.material = material
        self.boundary = boundary
        
    def CodeAsterExport(self, a=None):
        char_bound = ''
        for bound in self.boundary:
            char_bound += '_F(CHARGE = {},), '.format(bound.name)
        char = '''{} = MECA_STATIQUE(MODELE = {}, CHAM_MATER = {},
                EXCIT = ({}),); \n \n'''.format(self.name, self.modele.name, 
                self.material.name, char_bound)
        if a is None:
            return char
        else:
            a.write(char)
            
class ResultAnalysis:
    def __init__(self, simulation, typ, group_no):
        self.simulation = simulation
        self.typ = typ
        self.group_no = group_no
    
    def CodeAsterExport(self, a=None):
        char = '''{} = CALC_CHAMP(reuse = {},
                       RESULTAT = {},
                       FORCE = '{}',);
    
{} = POST_RELEVE_T(ACTION = _F(OPERATION = 'EXTRACTION',
                                 INTITULE = 'FORCE',
                                 RESULTAT = {},
                                 NOM_CHAM = '{}',
                                 GROUP_NO = '{}',
                                 RESULTANTE = ('DX', 'DY', 'DZ'),),);

IMPR_TABLE(TABLE = {}, UNITE = 4); \n \n'''.format(self.simulation.name,
        self.simulation.name, self.simulation.name, self.typ, self.simulation.name + '_1',
        self.simulation.name, self.typ, self.group_no.name, self.simulation.name + '_1')
        if a is None:
            return char
        else:
            a.write(char)
        

#IMPR_RESU(FORMAT='GMSH',
#          RESU=_F(MAILLAGE=MAIL,
#                  RESULTAT=RESU,
#                  NOM_CHAM='DEPL',),);
                
class CodeAster:
    def __init__(self, name, file_name, mesh, node_group, modele, material, 
                 boundary, simulation, analysis):
        
        self.name = name
        self.file_name = file_name
        self.mesh = mesh
        self.node_group = node_group
        self.modele = modele
        self.material = material
        self.boundary = boundary
        self.simulation = simulation
        self.analysis = analysis
        
    def GenerateCodeAsterFiles(self, path_code_aster):
        with open(self.file_name + '.export', 'w') as a:
            actual_path = os.getcwd()
            
            a.write('P actions make_etude\n')
            a.write('P aster_root {} \n'.format(path_code_aster))
            a.write('P consbtc oui\n')
            a.write('P corefilesize unlimited\n')
            a.write('P cpresok RESNOOK\n')
            a.write('P debug nodebug\n')
            a.write('P facmtps 1\n')
            a.write('P lang en\n')
            a.write('P memjob 524288\n')
            a.write('P memory_limit 1024.0\n')
            a.write('P mode interactif\n')
            a.write('P mpi_nbcpu 1\n')
            a.write('P mpi_nbnoeud 1\n')
            a.write('P nbmaxnook 5\n')
            a.write('P protocol_copyfrom asrun.plugins.server.SCPServer\n')
            a.write('P protocol_copyto asrun.plugins.server.SCPServer\n')
            a.write('P protocol_exec asrun.plugins.server.SSHServer\n')
            a.write('P soumbtc oui\n')
            a.write('P time_limit 900.0\n')
            a.write('P tpsjob 16\n')
            a.write('P version stable\n')
            a.write('A args\n')
            a.write('A memjeveux 64.0\n')
            a.write('A tpmax 900.0\n')
            
            a.write('F msh {}/{} D 19\n'.format(actual_path, self.mesh.file_name_mesh))
            a.write('F comm {}/{}.comm D 1\n'.format(actual_path, self.file_name))
            a.write('F resu {}/{}.resu R 8\n'.format(actual_path, self.file_name))
            a.write('F pos {}/{}.pos R 37\n'.format(actual_path, self.file_name))
            a.write('F mess {}/{}.mess R 6\n'.format(actual_path, self.file_name))
            a.write('F pos {}/{}_R.pos R 4\n'.format(actual_path, self.file_name))
            
            a.close()
        
    def CodeAsterExport(self):
        
        with open(self.file_name + '.comm', 'w') as a:
            a.write('DEBUT() \n \n')
            self.mesh.CodeAsterExport(a)
            for ng in self.node_group:
                ng.CodeAsterExport(a)
            self.modele.CodeAsterExport(a)
            self.material.CodeAsterExport(a)
            for bd in self.boundary:
                bd.CodeAsterExport(a)
            self.simulation.CodeAsterExport(a)
            for analyze in self.analysis:
                analyze.CodeAsterExport(a)
                
    #        char = '''IMPR_RESU(FORMAT='GMSH',
    #          {}=_F(MAILLAGE={},
    #                  RESULTAT={},
    #                  NOM_CHAM='DEPL',),); \n \n'''.format(self.simulation.name,
    #          self.mesh.name, self.simulation.name)
    #        a.write(char)
            a.write('FIN(); \n')
            a.close()
        
    def Run(self, path_code_aster='/opt/aster', version_code_aster='13.4'):
        self.CodeAsterExport()
        self.GenerateCodeAsterFiles(path_code_aster)
        
        actual_path = os.getcwd()
        os.system('cp {}/{}/share/aster/config.txt {}/.'.format(
                path_code_aster, version_code_aster, actual_path))
        
        os.system('{}/bin/./as_run {}.export'.format(path_code_aster, self.file_name))
        
    def Read(self, analysis, axis):
        actual_path = os.getcwd()
        with open(actual_path + '/' + self.file_name + '_R.pos', 'r') as a:
            lines = a.readlines()
            a.close()
            for line in lines:
                sp = line.split()
                if sp[0] == 'INTITULE':
                    liste_var = sp
                elif sp[0] == 'FORCE':
                    if axis == 'x':
                        pos = liste_var.index('DX')
                    elif axis == 'y':
                        pos = liste_var.index('DY')
                    elif axis == 'y':
                        pos = liste_var.index('DZ')
                    return float(sp[pos])
        
        
        