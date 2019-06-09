#!/usr/bin/env python

import os
import shutil
import numpy as np

from forcebalance.molecule import Molecule

molecules_folder = 'processed_molecules/mol2'
results_folder = 'td_results'
out_folder = 'targets'

target_str = '''
$target
name {name}
type AbInitio_SMIRNOFF
mol2 input.mol2
pdb conf.pdb
coords scan.xyz
writelevel 2
attenuate
energy_denom 2.0
energy_upper 10.0
force_rms_override 100.0
openmm_platform Reference
remote 1
$end
'''

def read_gradxyz(filename):
    m = Molecule()
    return m.read_xyz(filename)['xyzs']

def make_fb_targets():
    result_mol_folders = [os.path.join(results_folder, f) for f in os.listdir(results_folder) if os.path.isdir(os.path.join(results_folder, f))]
    result_mol_folders.sort()
    print(f"\nLoading data from {len(result_mol_folders)} result folders under {results_folder}")
    # output folder
    if os.path.exists(out_folder):
        shutil.rmtree(out_folder)
    os.mkdir(out_folder)
    target_names = []
    for mol_folder in result_mol_folders:
        mol_name = os.path.basename(mol_folder)
        # the name of the molecules should be consistent with the mol_folder
        mol_file = os.path.join(molecules_folder, mol_name + '.mol2')
        molecule = Molecule(mol_file)
        # find all torsion data
        finished_scans = []
        for f in os.listdir(mol_folder):
            name, ext = os.path.splitext(f)
            if ext == '.xyz':
                finished_scans.append(name)
        if len(finished_scans) == 0:
            print(f'No finished scans found in {mol_folder}')
            continue
        # output target name
        target_name = 'td_' + mol_name
        target_names.append(target_name)
        # make target folder
        this_target_folder = os.path.join(out_folder, target_name)
        os.mkdir(this_target_folder)
        # read data from each finished scans
        target_mol = Molecule()
        target_mol.elem = molecule.elem
        target_mol.xyzs = []
        target_mol.qm_energies = []
        target_mol.qm_grads = []
        for f in finished_scans:
            xyz_file = os.path.join(mol_folder, f+'.xyz')
            m = Molecule(xyz_file)
            target_mol.xyzs += m.xyzs
            # read energy from comment line
            energies = [float(comm.split()[-1]) for comm in m.comms]
            target_mol.qm_energies += energies
            # read gradient
            grad_file = os.path.join(mol_folder, f+'.gradxyz')
            grads = read_gradxyz(grad_file)
            target_mol.qm_grads += grads
        # write qdata.txt
        target_mol.write(os.path.join(this_target_folder, 'qdata.txt'))
        # write scan.xyz
        target_mol.write(os.path.join(this_target_folder, 'scan.xyz'))
        # write pdb
        molecule.write(os.path.join(this_target_folder, 'conf.pdb'))
        # copy mol2 file
        shutil.copyfile(mol_file, os.path.join(this_target_folder, 'input.mol2'))
        # write a note
        with open(os.path.join(this_target_folder, 'notes.txt'), 'w') as fnote:
            fnote.write("Notes: This target is made by make_fb_targets.py, using data from\n")
            fnote.write(mol_file + '\n')
            for f in finished_scans:
                xyz_file = os.path.join(mol_folder, f+'.xyz')
                grad_file = os.path.join(mol_folder, f+'.gradxyz')
                fnote.write(xyz_file + '\n')
                fnote.write(grad_file + '\n')
    # write a target.in file for use in ForceBalance input
    with open(os.path.join(out_folder, 'targets.in'), 'w') as fout:
        for tname in target_names:
            fout.write(target_str.format(name=tname) + '\n')
    print(f"Targets generation finished!")
    print(f"You can copy contents in {os.path.join(out_folder, 'targets.in')} to your ForceBalance input file.")

def main():
    if not os.path.isdir(molecules_folder):
        print(f"{molecules_folder} folder not found")
        return
    if not os.path.isdir(results_folder):
        print(f"{results_folder} folder not found")
        return
    make_fb_targets()

if __name__ == '__main__':
    main()
