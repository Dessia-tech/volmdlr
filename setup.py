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
    with open('README.md') as f:
        return f.read()
    
def version_from_git_describe(version):
    if version[0]=='v':
            version = version[1:]
            
    # PEP 440 compatibility
    number_commits_ahead = 0
    if '-' in version:
        version, number_commits_ahead, commit_hash = version.split('-')
        number_commits_ahead = int(number_commits_ahead)
        
    # print('number_commits_ahead', number_commits_ahead)
    
    split_versions = version.split('.')
    if 'post' in split_versions[-1]:
        suffix = split_versions[-1]
        split_versions = split_versions[:-1]
    else:
        suffix = None
    
    for pre_release_segment in ['a', 'b', 'rc']:
        if pre_release_segment in split_versions[-1]:
            if number_commits_ahead > 0:
                split_versions[-1] = str(split_versions[-1].split(pre_release_segment)[0])
                if len(split_versions) == 2:
                    split_versions.append('0')
                if len(split_versions) == 1:
                    split_versions.extend(['0', '0'])
        
                split_versions[-1] = str(int(split_versions[-1])+1)
                future_version = '.'.join(split_versions)
                return '{}.dev{}'.format(future_version, number_commits_ahead)
            else:
                return '.'.join(split_versions)
            
    if number_commits_ahead > 0:
        if len(split_versions) == 2:
            split_versions.append('0')
        if len(split_versions) == 1:
            split_versions.extend(['0', '0'])
        split_versions[-1] = str(int(split_versions[-1])+1)
        split_versions = '.'.join(split_versions)
        return '{}.dev{}'.format(split_versions, number_commits_ahead)
    else:
        if suffix is not None:
            split_versions.append(suffix)

        return '.'.join(split_versions)
    
# Just testing if get_version works well
assert version_from_git_describe('v0.1.7.post2') == '0.1.7.post2'
assert version_from_git_describe('v0.0.1-25-gaf0bf53') == '0.0.2.dev25'
assert version_from_git_describe('v0.1-15-zsdgaz') == '0.1.1.dev15'
assert version_from_git_describe('v1') == '1'
assert version_from_git_describe('v1-3-aqsfjbo') == '1.0.1.dev3'
    
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
        
        return version_from_git_describe(version)
    else:
        # Extract the version from the PKG-INFO file.
        with open(join(d, 'PKG-INFO')) as f:
            version = version_re.search(f.read()).group(1)
            
    # print('version', version)
    return version



setup(name='volmdlr',
      version = get_version(),
#      setup_requires=['setuptools_scm'],
      description=' A volume modeler computation-oriented. Include rendering bindings.',
      long_description=readme(),
      keywords='volume, modeler, CAD',
      url='https://github.com/Dessia-tech/volmdlr',
      author='DessiA Technologies',
      author_email='root@dessia.tech',
      license='Creative Commons Attribution-Share Alike license',
      packages=['volmdlr', 'volmdlr.models', 'volmdlr.models.workflows'],#,'volmdlr.primitives2D','volmdlr.primitives3D','volmdlr.geometry'],
      package_dir={},
      include_package_data = True,
      install_requires=['packaging',
                        'dessia_common>=0.3.9',
                        'numpy',
                        'matplotlib',
                        'Cython',
                        'scipy',
                        'geomdl',
                        'jsonschema',
                        'networkx',
                        'triangle',
                        'plot_data>=0.5.4'],
      classifiers=['Topic :: Scientific/Engineering','Development Status :: 3 - Alpha'],
      ext_modules = cythonize(["volmdlr/core_compiled.pyx"]))
