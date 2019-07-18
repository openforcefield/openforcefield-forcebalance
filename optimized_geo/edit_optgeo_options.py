#!/usr/bin/env python
# coding: utf-8

import os
folders = [f for f in os.listdir('.') if os.path.isdir(f)]
optfiles = [os.path.join(fol, 'optgeo_options.txt') for fol in folders]
optfiles = [f for f in optfiles if os.path.isfile(f)]

for f in optfiles:
    print(f"editing {f}")
    with open(f) as infile:
        lines = infile.readlines()
    for i, line in enumerate(lines):
        if line.startswith('dihedral_denom'):
            lines[i] = 'dihedral_denom 0.0\n'
    with open(f, 'w') as outfile:
        outfile.write(''.join(lines))
        
