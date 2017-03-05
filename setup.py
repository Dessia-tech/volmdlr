# -*- coding: utf-8 -*-
"""
Setup install script for volmdlr

"""

from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

import volmdlr
setup(name='volmdlr',
      version=volmdlr.__version__,#
      description=' A volume modeler computation-oriented. Include rendering bindings. ',
      long_description=readme(),
      keywords='volume, modeler',
      url='https://github.com/masfaraud/volmdlr',
      author='Steven Masfaraud',
      author_email='steven@masfaraud.fr',
      license='Creative Commons Attribution-Share Alike license',
      packages=['volmdlr'],#,'volmdlr.primitives2D','volmdlr.primitives3D','volmdlr.geometry'],
      package_dir={},
      install_requires=['numpy','matplotlib'],
      classifiers=['Topic :: Scientific/Engineering','Development Status :: 3 - Alpha'])
