"""
Sort the molecules by the change of objective from abinitio targets
"""
import os
import numpy as np
from sort_molecule_abinitio_obj import read_abinitio_notes, read_abinitio_EnergyCompare


def sort_print_obj_change_table(orig_molecule_obj_dict, new_molecule_obj_dict, mol_info):
    molecule_obj_change_dict = {m: new_molecule_obj_dict[m]-orig_molecule_obj_dict[m] for m in orig_molecule_obj_dict}
    sorted_molecules = sorted(orig_molecule_obj_dict, key=lambda x: molecule_obj_change_dict[x])
    print("Results")
    print(f"idx       {'folder':^70s} {'SMILES':^50s} objective")
    print('-'*120)
    for i, folder in enumerate(sorted_molecules):
        foldernm = 'targets/' + folder
        obj = molecule_obj_change_dict[folder]
        smiles = mol_info[folder]['smiles']
        print(f'{i:<7}   {foldernm:70s} {smiles:50s} {obj:9.5f}')


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--targets_folder', default='targets')
    parser.add_argument('-orig_tmp', '--orig_tmp_folder')
    parser.add_argument('-new_tmp', '--new_tmp_folder')
    args = parser.parse_args()

    mol_info = read_abinitio_notes(args.targets_folder)
    orig_molecule_obj_dict = read_abinitio_EnergyCompare(args.orig_tmp_folder)
    new_molecule_obj_dict = read_abinitio_EnergyCompare(args.new_tmp_folder)

    sort_print_obj_change_table(orig_molecule_obj_dict, new_molecule_obj_dict, mol_info)

if __name__ == "__main__":
    main()
