#!/usr/bin/env python

import os
import shutil
import subprocess
from geometric.molecule import Molecule
import numpy as np
import matplotlib.pyplot as plt

def read_scan_xyz(fnm):
    m = Molecule(fnm)
    grid_geo_data = {}
    grid_energies = {}
    for geo, comms in zip(m.xyzs, m.comms):
        dihedral_str = comms.split()[1][1:-1]
        grid_id = tuple(int(d) for d in dihedral_str.split(',') if d != '')
        grid_geo_data[grid_id] = geo
        grid_energies[grid_id] = float(comms.split()[-1])
    return grid_geo_data, grid_energies

def read_dihedralstxt(fnm):
    dihedral_list = []
    with open(fnm) as f:
        for line in f:
            line = line.split('#', maxsplit=1)[0].strip()
            if line:
                ls = line.split()
                if len(ls) >= 4:
                    dihedral_list.append([int(s) for s in ls[:4]])
    return dihedral_list

def write_constraints_txt(dihedral_list, grid_id):
    assert len(dihedral_list) == len(grid_id)
    with open('constraints.txt', 'w') as outfile:
        for dihedral, value in zip(dihedral_list, grid_id):
            outfile.write('$set\n')
            outfile.write('dihedral ' + ' '.join(map(str, dihedral)) + ' ' + str(value) + '\n')

tmpdir = 'constrain_opts_data'

def recompute_relax_torsion_profile(scanxyz, pdb_fnm, dihedralstxt, sysxml_fnm):
    grid_geo_data, _ = read_scan_xyz(scanxyz)
    dihedral_list = read_dihedralstxt(dihedralstxt)
    sysxml_path = os.path.abspath(sysxml_fnm)
    m = Molecule(pdb_fnm)
    # create a new tmp dir for running calculations
    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)
    os.mkdir(tmpdir)
    os.chdir(tmpdir)
    # loop over each geometry in each grid
    for grid_id, geo in grid_geo_data.items():
        folder = 'gid_' + '_'.join(f'{d:+03d}' for d in grid_id)
        print(f'Running constrained optimization in {folder}')
        os.mkdir(folder)
        os.chdir(folder)
        # copy files
        m.xyzs = [geo]
        m.write('frame.pdb')
        shutil.copy(sysxml_path, 'openmm_system.xml')
        write_constraints_txt(dihedral_list, grid_id)
        # run geometric
        command = "geometric-optimize openmm_system.xml constraints.txt --openmm --pdb frame.pdb --qccnv --reset --epsilon 0.0 --enforce 0.1 --qdata"
        subprocess.run(command, shell=True, check=True)
        os.chdir('..')
    os.chdir('..')

def plot_torsion_profile(scanxyz, plotfnm='mm_relaxed_energies.pdf'):
    _, grid_energies = read_scan_xyz(scanxyz)
    # data holder
    qm_energies = []
    mm_orig_energies = []
    mm_relaxed_energies = []
    # read mm energies from tmp folder
    os.chdir(tmpdir)
    sorted_grid_ids = sorted(grid_energies.keys())
    for grid_id in sorted_grid_ids:
        qm_energies.append(grid_energies[grid_id])
        # folder name consistent with recompute_relax_torsion_profile()
        folder = 'gid_' + '_'.join(f'{d:+03d}' for d in grid_id)
        m = Molecule(os.path.join(folder, 'qdata.txt'))
        initial_mm_energy = m.qm_energies[0]
        final_mm_energy = m.qm_energies[-1]
        mm_orig_energies.append(initial_mm_energy)
        mm_relaxed_energies.append(final_mm_energy)
    os.chdir('..')
    # compute relative energies in kcal/mol
    qm_energies = np.array(qm_energies)
    qm_energies = (qm_energies - qm_energies.min()) * 627.509
    mm_orig_energies = np.array(mm_orig_energies)
    mm_orig_energies = (mm_orig_energies - mm_orig_energies.min()) * 627.509
    mm_relaxed_energies = np.array(mm_relaxed_energies)
    mm_relaxed_energies = (mm_relaxed_energies - mm_relaxed_energies.min()) * 627.509
    # we can compute relaxed energies relative to the minimum of original mm for better understanding
    # mm_relaxed_energies = (mm_relaxed_energies - mm_orig_energies.min()) * 627.509
    # plot
    plt.Figure()
    plt.plot(sorted_grid_ids, qm_energies, label='QM Relative Energies')
    plt.plot(sorted_grid_ids, mm_orig_energies, label='MM Original Energies')
    plt.plot(sorted_grid_ids, mm_relaxed_energies, label='MM Relaxed Energies')
    plt.legend()
    plt.xlabel('Torsion Angle (degree)')
    plt.ylabel('Relative Eneriges (kcal/mol)')
    plt.savefig(plotfnm)
    print(f'Plot saved as {plotfnm}')



def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('scanxyz')
    parser.add_argument('dihedralstxt')
    parser.add_argument('pdb')
    parser.add_argument('sysxml')
    parser.add_argument('-p', '--plot_only', action='store_true', help='skip compute, just plot with existing data')
    args = parser.parse_args()

    if not args.plot_only:
        recompute_relax_torsion_profile(args.scanxyz, args.pdb, args.dihedralstxt, args.sysxml)
    plot_torsion_profile(args.scanxyz)

if __name__ == "__main__":
    main()