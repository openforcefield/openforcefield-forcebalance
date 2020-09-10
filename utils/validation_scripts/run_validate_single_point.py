#!/usr/bin/env python

import os
import pickle
import time

from simtk import openmm, unit

from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.topology import Molecule as OffMolecule
from openforcefield.topology import Topology as OffTopology


def eval_single_point_energies(forcefield, partial_charges, mol2_path, pdb_path, positions_list):
    """ Evaluate single-point energies of a molecule on one or more positions """
    off_mol = OffMolecule.from_file(mol2_path)
    # load partial charges
    off_mol.partial_charges = partial_charges
    off_top = OffTopology.from_molecules([off_mol])
    system = forcefield.create_openmm_system(off_top, charge_from_molecules=[off_mol])
    # create openmm simulation
    omm_topology = openmm.app.PDBFile(pdb_path).topology
    integrator = openmm.LangevinIntegrator(273.15, 1.0, 1.0)
    simulation = openmm.app.Simulation(omm_topology, system, integrator)
    # evaluate energy of each position
    result_energy_list = []
    for position in positions_list:
        simulation.context.setPositions(position)
        e = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        result_energy_list.append(e)
    return result_energy_list

def run_validate_single_point(ffxml, data_folder):
    forcefield = ForceField(ffxml, allow_cosmetic_attributes=True)
    ff_name = os.path.basename(ffxml)
    # read molecules from the data folder
    bench_mol_folders = sorted(os.listdir(data_folder))
    print(f'\n*** Single point benchmark on QM and MM minimized geometries ***')
    print(f"- Number of molecules: {len(bench_mol_folders)}")
    print(f"- Test forcefield : {ff_name}")
    print(f"- Start Time      : {time.ctime()}")
    print(f"- Unit            : kcal/mol \n")
    print(f'{"":5s} | {"":15s} | {"QM Geometry":^47s} | {"MM Geometry":^47s}')
    print(f'{"#":>5s} | {"Molecule ID":15s} | {"E_release":>15s} {"E_test":>15s} {"Diff":>15s} | {"E_release":>15s} {"E_test":>15s} {"Diff":>15s}')
    for i, mol_fol in enumerate(bench_mol_folders):
        folder_path = os.path.join(data_folder, mol_fol)
        # read bench data from folder
        with open(os.path.join(folder_path, 'bench_data.pickle'), 'rb') as pfile:
            bench_data = pickle.load(pfile)
        mol2_path = os.path.join(folder_path, bench_data['mol2_fnm'])
        pdb_path = os.path.join(folder_path, bench_data['pdb_fnm'])
        positions_list = [bench_data['geo1_qm'], bench_data['geo2_mm']]
        ref_energies = [bench_data['energy_geo1'], bench_data['energy_geo2']]
        eval_energies = eval_single_point_energies(forcefield, bench_data['partial_charges'], mol2_path, pdb_path, positions_list)
        e_ref_geo1 = ref_energies[0].value_in_unit(unit.kilocalorie_per_mole)
        e_ref_geo2 = ref_energies[1].value_in_unit(unit.kilocalorie_per_mole)
        e_test_geo1 = eval_energies[0].value_in_unit(unit.kilocalorie_per_mole)
        e_test_geo2 = eval_energies[1].value_in_unit(unit.kilocalorie_per_mole)
        print(f'{i:5d} | {mol_fol:15s} | {e_ref_geo1:15.3f} {e_test_geo1:15.3f} {e_test_geo1-e_ref_geo1:15.3f} | {e_ref_geo2:15.3f} {e_test_geo2:15.3f} {e_test_geo2-e_ref_geo2:15.3f}')


def main():
    import argparse
    parser = argparse.ArgumentParser('Prepare a validate folder containing data for single-point validations, from optgeo targets')
    parser.add_argument('ffxml', help='Forcefield file to benchmark')
    parser.add_argument('-d', '--data_folder', default='mol_data', help='Path to the data folder')
    args = parser.parse_args()

    run_validate_single_point(args.ffxml, args.data_folder)

if __name__ == "__main__":
    main()