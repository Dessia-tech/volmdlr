#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 22:50:39 2022

@author: steven
"""

# from radon.complexity import cc
import os
from pylama.main import check_paths, parse_options

MAX_RADON_SCORE = 20
RADON_HELP_MESSAGE = 'Use radon (https://radon.readthedocs.io/en/latest/commandline.html) to get some help on this code quality error.'

# Use and/or modify 0 or more of the options defined as keys in the variable my_redefined_options below.
# To use defaults for any option, remove that key completely.
my_redefined_options = {
    'linters': ['radon'],
    'skip': 'build/*,scripts',
    'verbose': 1,
}
# relative path of the directory in which pylama should check
# vm_path = '.'
vm_path = os.getcwd()
print(vm_path)
print(os.listdir(vm_path))

options = parse_options([vm_path], **my_redefined_options)
print('options', options)
errors = list(check_paths('volmdlr', options=options))

assert len(errors)

radon_raised = False
for error in errors:
    # print('{}: {}'.format(error.filename, error.message))
    if error.number == 'R901':
        radon_score = int(error.message.split(' ')[-1])
        if radon_score > MAX_RADON_SCORE:
            radon_raised = True
            print('Code quality error on {} @ line {}: {} (max {}): {}'.format(
                error.filename, error.lnum, radon_score, MAX_RADON_SCORE, error.message))

# print('radon raised', radon_raised)
if radon_raised:
    print(RADON_HELP_MESSAGE)
    raise RuntimeError('Radon error detected\n')
else:
    print('No radon error raised')