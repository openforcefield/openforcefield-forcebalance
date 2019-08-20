#!/usr/bin/env python

import os
import shutil
import subprocess
import pickle

from simtk import openmm

from forcebalance.molecule import Molecule

from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.topology import Molecule as OffMolecule
from openforcefield.topology import Topology as OffTopology

this_file_folder = os.path.dirname(os.path.realpath(__file__))

def read_molecules_from_optgeo_target(target_fol):
    # read all molecules from an optgeo target folder
    res = {}
    with open(os.path.join(target_fol, 'notes.txt')) as fnote:
        for line in fnote:
            ls = line.split()
            if len(ls) == 7 and ls[-1] == 'SUCCESS':
                assert ls[3] == 'molecule_id'
                mol_id = ls[4]
                fnm_pre = ls[0]
                # read content of mol2 file
                mol2_fnm = fnm_pre + '.mol2'
                with open(os.path.join(target_fol, mol2_fnm)) as fmol2:
                    mol2_str = fmol2.read()
                # read content of pdb file
                pdb_fnm = fnm_pre + '.pdb'
                with open(os.path.join(target_fol, pdb_fnm)) as f:
                    pdb_str = f.read()
                # read QM geometry from xyz file
                xyz_fnm = fnm_pre + '.xyz'
                m = Molecule(os.path.join(target_fol, xyz_fnm))
                data = {
                    'mol2_str': mol2_str,
                    'pdb_str': pdb_str,
                    'qm_geo': m.xyzs[0]
                }
                res[mol_id] = data
    return res

def gen_bench_data(forcefield, data):
    # write mol2 file
    mol2_fnm = 'input.mol2'
    with open(mol2_fnm, 'w') as f:
        f.write(data['mol2_str'])
    # write pdb file
    pdb_fnm = 'conf.pdb'
    with open(pdb_fnm, 'w') as f:
        f.write(data['pdb_str'])
    off_mol = OffMolecule.from_file(mol2_fnm)
    # compute the partial charges explicitly
    off_mol.generate_conformers()
    off_mol.compute_partial_charges_am1bcc()
    # save partial charges
    partial_charges = off_mol.partial_charges
    top = OffTopology.from_molecules([off_mol])
    system = forcefield.create_openmm_system(top, charge_from_molecules=[off_mol])
    # create simulation
    omm_topology = openmm.app.PDBFile(pdb_fnm).topology
    integrator = openmm.LangevinIntegrator(273.15, 1.0, 1.0)
    simulation = openmm.app.Simulation(omm_topology, system, integrator)
    # qm position and energy
    qm_pos = data['qm_geo'] * openmm.unit.angstrom
    simulation.context.setPositions(qm_pos)
    e_geo1 = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    # mm minimized position and energy
    for _ in range(10):
        # run the minimization several times to make it more stable
        simulation.minimizeEnergy(tolerance=1e-8, maxIterations=10000)
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    mm_pos = state.getPositions(asNumpy=True)
    e_geo2 = state.getPotentialEnergy()
    bench_data = {
        'mol2_fnm': mol2_fnm,
        'pdb_fnm': pdb_fnm,
        'partial_charges': partial_charges,
        'geo1_qm': qm_pos,
        'geo2_mm': mm_pos,
        'energy_geo1': e_geo1,
        'energy_geo2': e_geo2,
    }
    return bench_data

def gen_validate_single_point_folder(ffxml, targets_folder, targets_prefix):
    forcefield = ForceField(ffxml, allow_cosmetic_attributes=True)
    # collect molecule QM data from target folders
    mol_data = {}
    for target_fol in os.listdir(targets_folder):
        if target_fol.startswith(targets_prefix):
            path = os.path.join(targets_folder, target_fol)
            data = read_molecules_from_optgeo_target(path)
            print(f'Read {len(data)} molecule QM data from {path}')
            mol_data.update(data)
    print(f"Total {len(mol_data)} molecules found from target folders")
    # create a new folder
    fol_nm = 'validate_single_point'
    if os.path.exists(fol_nm):
        shutil.rmtree(fol_nm)
    os.mkdir(fol_nm)
    os.chdir(fol_nm)
    # create a "mol_data" folder to put all molecule data
    os.mkdir('mol_data')
    os.chdir('mol_data')
    # for each molecule, create a new subfolder and place the relative files in there
    print(f'Generating benchmark data by MM energy evaluation and minimizations')
    i = 0
    for mol_id, data in mol_data.items():
        os.mkdir(mol_id)
        os.chdir(mol_id)
        # generate benchmark data for a single molecule
        bench_data = gen_bench_data(forcefield, data)
        # save the bench_data as pickle file
        with open('bench_data.pickle', 'wb') as pfile:
            pickle.dump(bench_data, pfile)
        os.chdir('..')
        i += 1
        if i % 10 == 0:
            print(f'processed {i:4d} molecules')
    os.chdir('..')
    # copy and paste the "run_validate_single_point.py" script to here
    shutil.copyfile(os.path.join(this_file_folder, 'run_validate_single_point.py'), 'run_validate_single_point.py')
    os.chdir('..')

def main():
    import argparse
    parser = argparse.ArgumentParser('Prepare a validate folder containing data for single-point validations, from optgeo targets')
    parser.add_argument('-x', '--ffxml', help='Forcefield file to use')
    parser.add_argument('-t', '--targets_folder', help='Path to the FB targets folder')
    parser.add_argument('-p', '--targets_prefix', default='optgeo_SMIRNOFF_Coverage_Set', help='Prefix of targets for selection')
    args = parser.parse_args()

    gen_validate_single_point_folder(args.ffxml, args.targets_folder, args.targets_prefix)

if __name__ == "__main__":
    main()