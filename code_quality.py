#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 22:50:39 2022

@author: steven
"""

# from radon.complexity import cc

from pylama.main import check_paths, parse_options

# Use and/or modify 0 or more of the options defined as keys in the variable my_redefined_options below.
# To use defaults for any option, remove that key completely.
my_redefined_options = {
    'linters': ['radon'],
}
# relative path of the directory in which pylama should check
my_path = 'volmdlr'

options = parse_options([my_path], **my_redefined_options)
for error in check_paths(my_path, options, rootdir='.'):
    print(error)

