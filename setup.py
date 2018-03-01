# -*- coding: utf-8 -*-
"""
Setup install script for volmdlr

"""

from setuptools import setup
#from distutils.core import setup
from Cython.Build import cythonize


def readme():
    with open('README.rst') as f:
        return f.read()
    
def version_scheme(version):
    return '.'.join([str(i) for i in version.tag._key[1]])

def local_scheme(version):
    return ''

setup(name='volmdlr',
      use_scm_version={'version_scheme':version_scheme,'local_scheme':local_scheme},
      setup_requires=['setuptools_scm'],
      description=' A volume modeler computation-oriented. Include rendering bindings. ',
      long_description=readme(),
      keywords='volume, modeler',
      url='https://github.com/masfaraud/volmdlr',
      author='Steven Masfaraud',
      author_email='steven@masfaraud.fr',
      license='Creative Commons Attribution-Share Alike license',
      packages=['volmdlr'],#,'volmdlr.primitives2D','volmdlr.primitives3D','volmdlr.geometry'],
      package_dir={},
      install_requires=['numpy','matplotlib','Cython'],
      classifiers=['Topic :: Scientific/Engineering','Development Status :: 3 - Alpha'],
      ext_modules = cythonize("volmdlr/vmcy.pyx"))
