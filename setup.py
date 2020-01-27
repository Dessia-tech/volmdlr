# -*- coding: utf-8 -*-
"""
Setup install script for volmdlr

"""

from setuptools import setup
#from distutils.core import setup
from Cython.Build import cythonize

from os.path import dirname, isdir, join
import re
from subprocess import CalledProcessError, check_output

tag_re = re.compile(r'\btag: %s([0-9][^,]*)\b')
version_re = re.compile('^Version: (.+)$', re.M)

def readme():
    with open('README.rst') as f:
        return f.read()
    
def get_version():
    # Return the version if it has been injected into the file by git-archive
    version = tag_re.search('$Format:%D$')
    if version:
        return version.group(1)

    d = dirname(__file__)
    
    if isdir(join(d, '.git')):
        cmd = 'git describe --tags'
        try:
            version = check_output(cmd.split()).decode().strip()[:]
            
        except CalledProcessError:
            raise RuntimeError('Unable to get version number from git tags')
        if version[0]=='v':
            version = version[1:]
        # PEP 440 compatibility
        number_commits_ahead = 0
        if '-' in version:
            version, number_commits_ahead, commit_hash = version.split('-')
            number_commits_ahead = int(number_commits_ahead)

        future_version = version.split('.')

        for suffix in ['a', 'b', 'rc', 'post']:
            if suffix in future_version[-1]:
                if number_commits_ahead > 0:
                    future_version[-1] = str(future_version[-1].split(suffix)[0])
                    future_version[-1] = str(int(future_version[-1])+1)
                    future_version = '.'.join(future_version)
                    return '{}.dev{}'.format(future_version, number_commits_ahead)

                else:
                    return '.'.join(future_version)

        future_version[-1] = str(int(future_version[-1])+1)
        future_version = '.'.join(future_version)
        return '{}.dev{}'.format(future_version, number_commits_ahead)
        

    else:
        # Extract the version from the PKG-INFO file.
        with open(join(d, 'PKG-INFO')) as f:
            version = version_re.search(f.read()).group(1)
            
    print('version', version)
    return version


setup(name='volmdlr',
      version = get_version(),
#      setup_requires=['setuptools_scm'],
      description=' A volume modeler computation-oriented. Include rendering bindings. ',
      long_description=readme(),
      keywords='volume, modeler',
      url='https://github.com/Dessia-tech/volmdlr',
      author='DessiA Technologies',
      author_email='root@dessia.tech',
      license='Creative Commons Attribution-Share Alike license',
      packages=['volmdlr'],#,'volmdlr.primitives2D','volmdlr.primitives3D','volmdlr.geometry'],
      package_dir={},
      include_package_data = True,
      install_requires=['dessia_common', 'numpy', 'matplotlib', 'Cython', 'scipy', 'jinja2',
                        'geomdl', 'jsonschema', 'networkx', 'triangle'],
      classifiers=['Topic :: Scientific/Engineering','Development Status :: 3 - Alpha'],
      ext_modules = cythonize("volmdlr/core_compiled.pyx"))
