#!/usr/bin/env python

import os
import sys
import json

from forcebalance.molecule import Molecule
from submit_torsiondrives import TorsionSubmitter

def read_pick_csv(filename):
    lines = open(filename).readlines()
    res = []
    for line in lines[1:]:
        line = line.strip()
        ls = line.split(',')
        molecule_idx = int(ls[0])
        molecule_name = ls[1]
        dihedral = [int(i) for i in ls[2].split('-')]
        filename = f'input_molecules/{molecule_idx:03d}_{molecule_name}.mol2'
        res.append([filename, dihedral])
    return res



def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", help='Input csv file containing definition of dihedrals')
    parser.add_argument("-s", "--scan_config", help='File containing configuration of QM scans')
    args = parser.parse_args()

    print(' '.join(sys.argv))

    submitter = TorsionSubmitter(scan_conf_file=args.scan_config)

    picked_data = read_pick_csv(args.csv)

    for f, dihedral in picked_data:
        print(f"Submitting {f} {dihedral}")
        submitter.submit_molecule(f, dihedral_list=[dihedral], to_json=True)

    submitter.write_submitted_json("picked_torsion_options.json")


if __name__ == '__main__':
    main()
