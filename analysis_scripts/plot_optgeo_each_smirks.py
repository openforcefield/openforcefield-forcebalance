#!/usr/bin/env python
"""
This script is used to analyze ForceBalance fitting resutls for optgeo targets.
Author: Yudong Qiu
"""

import os
import sys
import shutil
import numpy as np
import pickle

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from forcebalance.molecule import Molecule

from simtk import openmm, unit

from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.topology import Molecule as OffMolecule
from openforcefield.topology import Topology as OffTopology

FFTYPE_MAP = {
    'Distance': 'Bonds',
    'Angle': 'Angles',
    'Dihedral': 'ProperTorsions',
    'Out-of-Plane': 'ImproperTorsions',
}

def read_aggregate_optgeo_data(tmp_folder, iter_folder, forcefield, targets_folder='targets', compare_iter0=True):
    """ Read optgeo target data from tmp folder
    tmp_folder: optimize.tmp
    iter_folder: iter_0040
    """
    optgeo_folders = [os.path.join(tmp_folder, f) for f in os.listdir(tmp_folder) if 'optgeo' in f and os.path.isdir(os.path.join(tmp_folder, f))]
    optgeo_folders.sort()
    print(f"Reading optgeo target rmsd_decomposition from {len(optgeo_folders)} folders")
    data_qm_v_mm = {t:{} for t in FFTYPE_MAP.values()}
    for i, fol in enumerate(optgeo_folders):
        print(f"{i}: {fol}")
        rmsd_fnm = os.path.join(fol, iter_folder, 'rmsd_decomposition.txt')
        print(f" - reading {rmsd_fnm}")
        data = read_rmsd_decomposition(rmsd_fnm)
        if compare_iter0:
            rmsd_fnm_iter0 = os.path.join(fol, 'iter_0000', 'rmsd_decomposition.txt')
            print(f" - reading {rmsd_fnm_iter0}")
            data_iter0 = read_rmsd_decomposition(rmsd_fnm_iter0)
        print(f" - found data for {len(data)} molecules")
        print(f" - parsing smirks for each molecule")
        tgt_folder = os.path.join(targets_folder, os.path.basename(fol))
        for mol_name in data.keys():
            print(f"   - {mol_name}")
            mol2_fnm = os.path.join(tgt_folder, mol_name+'.mol2')
            smirnoff_assignments = compute_smirnoff_assignments(forcefield, mol2_fnm)
            for fftype in data[mol_name]:
                for atom_indices, (qm, mm) in data[mol_name][fftype].items():
                    try:
                        sid, smirks = smirnoff_assignments[fftype][atom_indices]
                    except KeyError:
                        print(f"    Warning: {mol2_fnm} {fftype} {atom_indices} no match")
                        continue
                    data_qm_v_mm[fftype].setdefault(sid, [])
                    data_qm_v_mm[fftype][sid].append({
                        'mol2_fnm': mol2_fnm,
                        'atom_indices': atom_indices,
                        'smirks': smirks,
                        'id': sid,
                        'qm': qm,
                        'mm': mm,
                        'mm_iter0': data_iter0[mol_name][fftype][atom_indices][1] if compare_iter0 else None
                    })
    return data_qm_v_mm


def read_rmsd_decomposition(fnm):
    """
    {
        '001_CH4.mol2': {
            "Bonds": {
                (0,1): (1.5, 1.4),
                (0,2): (1.4, 1.3),
            },
            "Angles": {
                (0,1,2): (150, 148),
            },
            "ProperTorsions": {
                ...
            },
            ..
        }
    }
    """
    with open(fnm) as f:
        content = f.read()
    mol_blocks = content.split('\n\n')
    res = {}
    for block in mol_blocks:
        lines = block.split('\n')
        mol_name = None
        for line in lines:
            if not line: continue
            ls = line.split()
            if line.startswith('['):
                mol_name = ls[1]
                res[mol_name] = {t:{} for t in FFTYPE_MAP.values()}
            elif ls[0] in FFTYPE_MAP:
                fftype = FFTYPE_MAP[ls[0]]
                atom_idxs = tuple(int(a)-1 for a in ls[1].split('-'))
                # make ordering consistent between geometric and smirnoff
                if fftype == 'ImproperTorsions':
                    center, rest = atom_idxs[0], sorted(atom_idxs[1:])
                    atom_idxs = (rest[0], center, rest[1], rest[2])
                elif atom_idxs > atom_idxs[::-1]:
                    atom_idxs = atom_idxs[::-1]
                qm, mm, diff = float(ls[2]), float(ls[3]), float(ls[4])
                # use the minimum diff to replace mm
                mm = qm - diff
                # convert radian to degree
                # if fftype != 'Bonds':
                #     qm = qm * 180 / np.pi
                #     mm = mm * 180 / np.pi
                # if mol_name.startswith('350'):
                #     print(f"{mol_name} {fftype} {atom_idxs} {qm} {mm}")
                res[mol_name][fftype][atom_idxs] = (qm, mm)
    return res


def compute_smirnoff_assignments(forcefield, mol2_fnm):
    off_mol = OffMolecule.from_file(mol2_fnm)
    off_top = OffTopology.from_molecules(off_mol)
    label_result = forcefield.label_molecules(off_top)[0]
    res = {t:{} for t in FFTYPE_MAP.values()}
    # read bonds
    for atom_indices, bond_param in label_result['Bonds'].items():
        # equilibrium bond length in angstrom
        res['Bonds'][atom_indices] = (bond_param.id, bond_param.smirks)
    # read angles
    for atom_indices, angle_param in label_result['Angles'].items():
        # equilibrium angle in degree
        res['Angles'][atom_indices] = (angle_param.id, angle_param.smirks)
    # read dihiedrals
    for atom_indices, torsion_param in label_result['ProperTorsions'].items():
        res['ProperTorsions'][atom_indices] = (torsion_param.id, torsion_param.smirks)
    # read impropers
    for atom_indices, imp_torsion_param in label_result['ImproperTorsions'].items():
        res['ImproperTorsions'][atom_indices] = (imp_torsion_param.id, imp_torsion_param.smirks)
    return res

def generate_analysis_plots(data_qm_v_mm, orig_equilibrium, new_equilibrium, iter_folder):
    folder = 'optgeo_scatter_plots'
    if os.path.exists(folder):
        shutil.rmtree(folder)
    os.mkdir(folder)
    os.chdir(folder)
    # generate plots for each smirks pattern
    for fftype in data_qm_v_mm:
        print(f"{fftype}")
        os.mkdir(fftype)
        os.chdir(fftype)
        ff_data = data_qm_v_mm[fftype]
        for sid, dlist in ff_data.items():
            if len(dlist) > 0:
                qm_array = [d['qm'] for d in dlist]
                mm_array = [d['mm'] for d in dlist]
                mm_array_iter0 = [d['mm_iter0'] for d in dlist] if dlist[0]['mm_iter0'] != None else None
                smirks = dlist[0]['smirks']
                title = f"<{fftype}> {sid}: {smirks} [n={len(dlist)}]"
                plot_qm_mm_scatter(qm_array, mm_array, f"{sid}.pdf", scatter_label=iter_folder, mm_array_iter0=mm_array_iter0, title=title,
                    orig_equil=orig_equilibrium.get(sid), new_equil=new_equilibrium.get(sid))
        os.chdir('..')
    os.chdir('..')

def plot_qm_mm_scatter(qm_array, mm_array, fnm, scatter_label=None, mm_array_iter0=None, title="", orig_equil=None, new_equil=None):
    plt.Figure()
    if mm_array_iter0 != None:
        plt.scatter(mm_array_iter0, qm_array, marker='x', color='C1', alpha=0.5, label='iter_0')
    else:
        mm_array_iter0 = []
    plt.scatter(mm_array, qm_array, marker='x', alpha=0.5, color='C0', label=scatter_label)
    plt.xlabel('MM Value')
    plt.ylabel('QM Value')
    xmin, xmax = min(mm_array + mm_array_iter0), max(mm_array + mm_array_iter0)
    ymin, ymax = min(qm_array), max(qm_array)
    vmin, vmax = min(xmin, ymin), max(xmax, ymax)
    rng = vmax - vmin
    # eq value line
    if orig_equil != None and new_equil != None:
        vmin = min(vmin, orig_equil, new_equil)
        vmax = max(vmax, orig_equil, new_equil)
        plt.plot([orig_equil, orig_equil], [vmin, vmax], color='C1', alpha=0.5, label='orig equilibrium value')
        plt.plot([new_equil, new_equil], [vmin, vmax], color='C2', alpha=0.5, label='new equilibrium value')
    # reference diagnoal line
    plt.plot([vmin,vmax],[vmin,vmax], '--', color='black', alpha=0.5, label='reference QM=MM')
    plt.legend(framealpha=0.5)
    plt.title(title)
    # plt.tight_layout()
    plt.axis('equal')
    plt.axis([vmin-0.05*rng, vmax+0.05*rng, vmin-0.05*rng, vmax+0.05*rng])
    fig = plt.gcf()
    fig.set_size_inches(5,5)
    fig.savefig(fnm)
    # plt.savefig(fnm)
    plt.close()

def get_equilibrium_values(forcefield):
    equilibrium_values = {}
    for fftype in ['Bonds', 'Angles']:
        handler = forcefield.get_parameter_handler(fftype)
        for param in handler.parameters:
            if fftype == 'Bonds':
                v = param.length._value
            elif fftype == 'Angles':
                v = param.angle._value
            else:
                v = None
            equilibrium_values[param.id] = v
    return equilibrium_values


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-x', '--ffxml', default='forcefield/smirnoff99Frosst_experimental.offxml', help='Forcefield file to use')
    parser.add_argument('--new_xml', default='result/optimize/param_valence.offxml', help='New force field file after fitting')
    parser.add_argument('-i', '--iter', type=int, default=None, help='iteration number to read rmsd')
    parser.add_argument('-l', '--load_pickle', help='Load data directly from pickle file')

    args = parser.parse_args()

    if args.load_pickle:
        data_save = pickle.load(open(args.load_pickle, 'rb'))
        data_qm_v_mm = data_save['data_qm_v_mm']
        iter_folder = data_save['iter_folder']
    else:
        forcefield = ForceField(args.ffxml, allow_cosmetic_attributes=True)
        tmp_folder = 'optimize.tmp'
        if args.iter != None:
            iter_folder = f"iter_{args.iter:04d}"
        else:
            optgeo_folders = [os.path.join(tmp_folder, f) for f in os.listdir(tmp_folder) if 'optgeo' in f and os.path.isdir(os.path.join(tmp_folder, f))]
            iter_folder = max([f for f in os.listdir(optgeo_folders[0]) if f.startswith('iter_')])
        compare_iter0 = (iter_folder != 'iter_0000')
        data_qm_v_mm = read_aggregate_optgeo_data(tmp_folder, iter_folder, forcefield, compare_iter0=compare_iter0)
        # save the data on disk
        res_data_fnm = 'optgeo_analysis_data.p'
        with open(res_data_fnm, 'wb') as pfile:
            data_save = {
                'ffxml': args.ffxml,
                'tmp_folder': tmp_folder,
                'iter_folder': iter_folder,
                'data_qm_v_mm': data_qm_v_mm,
            }
            pickle.dump(data_save, pfile)

    # get eq values
    orig_forcefield = ForceField(args.ffxml, allow_cosmetic_attributes=True)
    new_forcefield = ForceField(args.new_xml, allow_cosmetic_attributes=True)
    orig_equilibrium = get_equilibrium_values(orig_forcefield)
    new_equilibrium = get_equilibrium_values(new_forcefield)

    # generate plots
    print("Generating plots")
    generate_analysis_plots(data_qm_v_mm, orig_equilibrium, new_equilibrium, iter_folder)

if __name__ == '__main__':
    main()