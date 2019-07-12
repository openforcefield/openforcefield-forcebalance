#!/usr/bin/env python

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

os.environ['OMP_NUM_THREADS'] = '1'

def build_context(offxml, molfile):
    """ Build an OpenMM Context from a offxml file and a molecule file """
    forcefield = ForceField(offxml, allow_cosmetic_attributes=True)
    molecule = OffMolecule.from_file(molfile)
    system = forcefield.create_openmm_system(molecule.to_topology())
    integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
    platform = openmm.Platform.getPlatformByName('Reference')
    context = openmm.Context(system, integrator, platform)
    return context

def evaluate_energies(context, trajectory):
    """ Evaluate energies of all frames in a trajectory using OpenMM context """
    energies = []
    for frame_geo in trajectory:
        context.setPositions(frame_geo * unit.angstrom)
        energy = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
        energies.append(energy)
    return np.array(energies)

def get_td_angles(fb_molecule):
    td_angles = []
    for comm in fb_molecule.comms:
        td_angles.append(int(comm.split()[1][1:-2]))
    return td_angles

def get_qm_energies(fb_molecule):
    """ Get qm energy information from fb_molecule, return relative energies in kcal/mol """
    energies = []
    for comm in fb_molecule.comms:
        energies.append(float(comm.rsplit(maxsplit=1)[-1]))
    energies = np.array(energies)
    energies -= energies.min()
    energies *= 627.509
    return energies

def plot_energies_data(energies_data_dict, filename):
    """ Generate an energy data plot fromt data dict """
    plt.Figure()
    x_axis = energies_data_dict.pop('td_angles')
    for dataname, datavalues in energies_data_dict.items():
        plt.plot(x_axis, datavalues, label=dataname)
    plt.legend()
    plt.title("Relative Energies Comparison")
    plt.ylabel("Relative Energies [ kcal/mol ]")
    plt.xlabel("Torsion Angle [ degree ]")
    plt.savefig(filename)
    plt.close()


def plot_torsion_energy_compare(offxml_list, td_results_folder, processed_molecules_folder):
    ff_names = [os.path.splitext(os.path.basename(f))[0] for f in offxml_list]
    target_molecule_folders = [os.path.join(td_results_folder, fol) for fol in os.listdir(td_results_folder) if os.path.isdir(os.path.join(td_results_folder, fol))]
    target_molecule_folders.sort()
    # create an empty output folder
    output_folder = 'td_energies_plots'
    if os.path.isdir(output_folder):
        shutil.rmtree(output_folder)
    os.mkdir(output_folder)
    for mol_folder in target_molecule_folders:
        print(f"{mol_folder}")
        # find mol file that contains topology
        mol_name = os.path.basename(mol_folder)
        mol_file = os.path.join(processed_molecules_folder, mol_name+'.sdf')
        assert os.path.isfile(mol_file), f'Molecule file {mol_file} not found'
        # build openmm contexts for each force field file
        contexts = [build_context(offxml, mol_file) for offxml in offxml_list]
        # find scan trajectories
        traj_files = [f for f in os.listdir(mol_folder) if os.path.splitext(f)[-1] == '.xyz']
        # create output subfolder
        out_mol_folder = os.path.join(output_folder, mol_name)
        os.mkdir(out_mol_folder)
        for f in traj_files:
            print(f"- {f}")
            # hold energy data for this trajectory
            energies_data_dict = {}
            # use ForceBalance.molecule to read the xyz file
            fb_mol = Molecule(os.path.join(mol_folder, f))
            # read torsion angles
            energies_data_dict['td_angles'] = get_td_angles(fb_mol)
            # read QM energies
            energies_data_dict['QM'] = eqm = get_qm_energies(fb_mol)
            # record the index of the minimum energy structure
            ground_idx = np.argmin(eqm)
            # evalute mm energies for each force field
            for context, ffname in zip(contexts, ff_names):
                mm_energies = evaluate_energies(context, fb_mol.xyzs)
                mm_energies -= mm_energies[ground_idx]
                energies_data_dict[ffname] = mm_energies
            # save the data on disk
            data_file_name = os.path.splitext(f)[0] + '_energies.pickle'
            with open(os.path.join(out_mol_folder, data_file_name), 'wb') as picklefile:
                pickle.dump(energies_data_dict, picklefile)
            # plt the data
            plot_file_name = os.path.splitext(f)[0] + '.pdf'
            plot_energies_data(energies_data_dict, os.path.join(out_mol_folder, plot_file_name))
    print(f"Evaluation finished, all data and plots saved in {output_folder}")

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('offxml', nargs='+', help='SMIRNOFF forcefield offxml files to evaluate energies')
    parser.add_argument('--td_results_folder', default='td_results', help='Folder name for the td_result data')
    parser.add_argument('--processed_molecules_folder', default='processed_molecules', help='Folder name for the processed molecules')
    args = parser.parse_args()

    plot_torsion_energy_compare(args.offxml, args.td_results_folder, args.processed_molecules_folder)

if __name__ == "__main__":
    main()