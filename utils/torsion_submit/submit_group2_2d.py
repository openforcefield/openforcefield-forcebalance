#!/usr/bin/env python

import os
import sys
import copy
import yaml
import json
import numpy as np


from process_molecules import read_sdf_to_fb_mol
from find_dihedrals import DihedralSelector
from torsion_submitter import TorsionSubmitter


def get_group2_2d_dihedrals(filename):
    m = read_sdf_to_fb_mol(filename)
    dihedral_selector = DihedralSelector(m)
    dihedral_pairs_list = dihedral_selector.find_dihedral_pairs(pattern='ring-a-ring')
    return dihedral_pairs_list

def submit_group2_2d(filenames, scan_conf_file, client_conf_file, to_json):
    submitter = TorsionSubmitter(scan_conf_file=scan_conf_file, client_conf_file=client_conf_file)
    for f in filenames:
        print(f"\n*** Submitting 2-D torsion scans for {f} ***")
        dihedral_pairs_list = get_group2_2d_dihedrals(f)
        submitter.submit_2d(f, dihedral_pairs_list)
    submitter.write_checkpoint()

def prepare_group2_2d_json(filenames, scan_conf_file):
    submitter = TorsionSubmitter(scan_conf_file=scan_conf_file)
    for f in filenames:
        print(f"\n*** Preparing 2-D torsion scans as JSON for {f} ***")
        dihedral_pairs_list = get_group2_2d_dihedrals(f)
        submitter.prepare_2d_json(f, dihedral_pairs_list)
    submitter.write_checkpoint()
    submitter.write_submitted_json("submit_torsion_options.json")


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("infiles", nargs='*', help='Input sdf file for a single molecule')
    parser.add_argument("-s", "--scan_config", help='File containing configuration of QM scans')
    parser.add_argument("-c", "--client_config", help='File containing configuration of QCFractal Client')
    parser.add_argument("-j", "--save_json", action="store_true", default=False, help='If specified, will not submit but save all submit job in a json.')
    args = parser.parse_args()

    print(' '.join(sys.argv))
    if args.save_json:
        prepare_group2_2d_json(args.infiles, args.scan_config)
    else:
        submit_group2_2d(args.infiles, args.scan_config, args.client_config)

if __name__ == '__main__':
    main()
