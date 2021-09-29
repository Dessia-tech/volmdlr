#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 14:35:47 2021

@author: steven
"""

import json

with open('coverage.json', 'r') as file:
    d = json.load(file)

assert d['totals']['percent_covered'] > 25.

print('total covered', d['totals']['percent_covered'], '%')

for file_name, data in d['files'].items():
    print(file_name, data['summary']['percent_covered'], '%')
    assert(data['summary']['percent_covered']) > 20.