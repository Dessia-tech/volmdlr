#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 14:07:37 2017

@author: steven
"""

class VolumeModel:
    def __init__(self,primitives):
        self.primitives=primitives
        
    def FreeCADScript(self,fcstd_filepath,path_lib_freecad='',py_filepath=''):
        """
        Generate python a FreeCAD definition of model 
        :param fcstd_filename: a filename without extension to give the name at the fcstd part written in python code
        :param py_filepath: python file to be exported. If none give, only return in string
        :type fcstd_filename:str
        :type py_filepath:str
        """
        s=''
        if path_lib_freecad!='':
            s+="import sys\nsys.path.append('"+path_lib_freecad+"')\n"

        s+="import FreeCAD as fc\nimport Part\n\ndoc_bv=fc.newDocument('bv')\n\n"
        
        for ip,primitive in enumerate(self.primitives):
            print(ip)
            s+=(primitive.FreeCADExport(ip)+'\n')
            s+=('shapeobj = doc_bv.addObject("Part::Feature","primitive'+str(ip)+': '+primitive.name+'")\n')
            s+=("shapeobj.Shape = primitive"+str(ip)+'\n\n')
            
        s+='doc_bv.recompute()\n'
        s+="doc_bv.saveAs('"+fcstd_filepath+".fcstd')"

        if py_filepath!='':
            with open(py_filepath,'w') as fichier:
                fichier.write(s)
        return s
            
    def FreeCADExport(self,python_path,fcstd_filepath,path_lib_freecad=''):
        """
        Export model to .fcstd FreeCAD standard
        
        :param python_path: path of python binded to freecad
             - on windows: 'something like C:\Program Files\FreeCAD X.XX\bin\python'
             - on linux: python (in general)
        :param filepath: path of fcstd file (without extension)
        :param path_lib_freecad: FreeCAD.so lib path (/usr/lib/freecad/lib in general)

        """
        s=self.FreeCADScript(fcstd_filepath,path_lib_freecad)
        print(s)
        import subprocess
        arg='-c\n'+s
        rep=subprocess.call([python_path,arg])
        print(rep)