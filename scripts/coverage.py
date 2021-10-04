#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 14:35:47 2021

@author: steven
"""

import json

untracked_modules = ['volmdlr/templates.py',
                     'volmdlr/code_aster.py',
                     'volmdlr/core_compiled.py',
                     'volmdlr/mesh.py',
                     'models/__init__.py',
                     'workflows/__init__.py',
                     'workflows/core.py',
                     'volmdlr/cloud.py']

with open('coverage.json', 'r') as file:
    d = json.load(file)


print('total covered', d['totals']['percent_covered'], '%')
assert d['totals']['percent_covered'] > 25.

for file_name, data in d['files'].items():
    print(file_name, data['summary']['percent_covered'], '%')
    print('/'.join(file_name.split('/')[-2:])) 
    if '/'.join(file_name.split('/')[-2:]) in untracked_modules:
        print(file_name, '-> in untrack list')
    else:
        assert(data['summary']['percent_covered']) > 20.
