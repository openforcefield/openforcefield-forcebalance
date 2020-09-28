#!/usr/bin/env python

import os
import shutil
import pickle
import numpy as np
import copy
import json

import cmiles
from openeye import oechem

from forcebalance.molecule import Molecule
from openforcefield.typing.engines.smirnoff import ForceField

from simtk.openmm import app, unit

import qcportal as ptl

client = ptl.FractalClient('https://api.qcarchive.molssi.org:443/')
ofs = oechem.oemolostream()


def download_hessian_data(dataset_name):
    """
    Download data from public server

    Parameters
    ----------
    dataset_name: str
        example: "OpenFF Optimization Set 1"

    Returns
    -------
    hessian_data: dict
        {
            'COC(O)OC-0': {
                'molecule': <Molecule(name='C6H8O' formula='C6H8O' hash='8c27507')>,
                'energy': -308.71238270980876,
                'gradient': np.array([[4.90e-05, 3.98e-05, 2.37e-05], ...]), # shape N x 3
                'hessian': np.array([[0.624, -0.030, 0.035, ...], ...]), # shape 3N x 3N
                'qcvars': {
                    '-D ENERGY': -0.02057017,
                    'B3LYP-D3(BJ) DISPERSION CORRECTION ENERGY': -0.02057017,
                    'CURRENT DIPOLE X': -0.019230813890716552,
                    'CURRENT DIPOLE Y': 0.29464541067283934,
                    'CURRENT DIPOLE Z': 0.5339113783843144,
                    'CURRENT ENERGY': -308.71238270980876,
                    'CURRENT REFERENCE ENERGY': -308.71238270980876,
                    'DFT FUNCTIONAL TOTAL ENERGY': -308.6918125398088,
                    'DFT TOTAL ENERGY': -308.71238270980876,
                    'DFT VV10 ENERGY': 0.0,
                    'DFT XC ENERGY': -35.58757022575535,
                    'DISPERSION CORRECTION ENERGY': -0.02057017,
                    'NUCLEAR REPULSION ENERGY': 288.5270197383624,
                    'ONE-ELECTRON ENERGY': -988.2068787424026,
                    'PCM POLARIZATION ENERGY': 0.0,
                    'SCF DIPOLE X': -0.019230813890716552,
                    'SCF DIPOLE Y': 0.29464541067283934,
                    'SCF DIPOLE Z': 0.5339113783843144,
                    'SCF ITERATION ENERGY': -308.71238270980876,
                    'SCF ITERATIONS': 12.0,
                    'TWO-ELECTRON ENERGY': 426.5756166899869,
                    'CURRENT GRADIENT': [...],
                    ...
                },
            },
            ...
        }
    """
    # load dataset from public qcfractal server
    # temporary: We don't have a "HessianDataset" yet, so many operations are "manual" on the general Dataset
    ds = client.get_collection("Dataset", dataset_name)
    # read specs
    spec_values = list(ds.data.history)[0]
    spec_dict = dict(zip(ds.data.history_keys, spec_values))
    assert spec_dict['driver'] == 'hessian'
    method = spec_dict['method']
    basis = spec_dict['basis']
    print(f"Specs for [ {dataset_name} ] loaded\n{spec_dict}")
    # download data for all molecules
    # dict_mol_id_entry_name = {rec.molecule_id: rec.name for rec in ds.data.records}
    dict_mol_id_entry_name = {row["molecule_id"]: row["name"] for index, row in ds.get_entries().iterrows()}
    print(f"Found total {len(dict_mol_id_entry_name)} molecule entries")
    all_mol_ids = list(dict_mol_id_entry_name.keys())
    entry_molecule_dict = {}
    for mol in client.query_molecules(all_mol_ids):
        entry_name = dict_mol_id_entry_name[mol.id]
        entry_molecule_dict[entry_name] = mol
    # query compute record
    hessian_data = {}
    print(f"Downloading hessian job data for {len(all_mol_ids)} molecules")
    for record in client.query_results(molecule=all_mol_ids, driver="hessian", method=method, basis=basis):
        if record.status != 'COMPLETE':
            continue
        mol_id = record.molecule
        entry_name = dict_mol_id_entry_name[mol_id]
        molecule = entry_molecule_dict[entry_name]
        noa = len(molecule.symbols)
        gradient = np.array(record.extras['qcvars']['CURRENT GRADIENT'], dtype=float).reshape(noa, 3)
        hessian = np.array(record.return_result, dtype=float).reshape(noa*3, noa*3)
        hessian_data[entry_name] = {
            'molecule': molecule,
            'energy': record.properties.return_energy,
            'gradient': gradient,
            'hessian': hessian,
            'qcvars': record.extras['qcvars'],
        }
    print(f'Downloaded hessian data for {len(hessian_data)} completed entries')
    # save as pickle file
    with open('hessian_data.pickle', 'wb') as pfile:
        pickle.dump(hessian_data, pfile)
    return hessian_data

def filter_hessian_data(hessian_data):
    """ Filter contents of the hessian data, pick the entry with lowest energy """
    res = {}
    lowest_energy = {}
    print("Start filtering hessian data by picking the entry with lowest energy")
    for entry_name, data in hessian_data.items():
        mol_name = entry_name.rsplit('-', maxsplit=1)[0]
        if mol_name not in lowest_energy:
            lowest_energy[mol_name] = data['energy']
            res[entry_name] = data
        else:
            current_energy = lowest_energy[mol_name]
            if data['energy'] < current_energy:
                lowest_energy[mol_name] = data['energy']
                res[entry_name] = data
    print(f"Filter hessian data complete, {len(res)} data entries left")
    return res

target_in_str = '''
$target
name {name}
type VIBRATION_SMIRNOFF
coords conf.pdb
mol2 input.mol2
weight 1.0
wavenumber_tol 200.0
remote 1
$end
'''

def make_vib_freq_target(dataset_name, hessian_data, test_ff=None):
    """
    Make a folder of ForceBalance targets from the torsiondrive data
    """
    target_name_prefix = 'vibfreq_' + dataset_name.replace(' ', '_')
    # create new targets folder
    if os.path.exists('targets'):
        shutil.rmtree('targets')
    os.mkdir('targets')
    os.chdir('targets')
    # write each entry as an individual target
    target_idx = 0
    n_targets = len(hessian_data)
    idx_fmt_string = get_int_fmt_string(n_targets)
    target_names = []
    print(f"Generating {n_targets} targets")
    for entry_index, data in hessian_data.items():
        # get formula for the molecule
        qcmol = data['molecule']
        # get mol_formula
        mol_formula = qcmol.get_molecular_formula()
        # create target folder
        target_idx_str = idx_fmt_string.format(target_idx)
        target_name = f"{target_name_prefix}_{target_idx_str}_{mol_formula}"
        print(f"{target_idx_str}: {target_name:70s} - ", end='')
        os.mkdir(target_name)
        os.chdir(target_name)
        # write a note
        with open('note.txt', 'w') as notefile:
            notefile.write(f'Target generated from dataset {dataset_name}, entry {entry_index}')
        # write input.mol2 file
        qcjson_mol = qcmol.dict(encoding='json')
        success = True
        try:
            oemol = cmiles.utils.load_molecule(qcjson_mol)
            ofs.open('input.mol2')
            oechem.OEWriteMolecule(ofs, oemol)
            ofs.close()
            # test if bonds changed
            if not check_connectivity('input.mol2'):
                success = False
                msg = "Bonds changed after rebuild"
            if success and test_ff != None:
                success, msg, molecule_labels = test_ff_mol2(test_ff, 'input.mol2')
            if not success:
                if not os.path.exists('../error_mol2s'):
                    os.mkdir('../error_mol2s')
                shutil.move(f'input.mol2', f'../error_mol2s/{target_name}.mol2')
                with open(f'../error_mol2s/{target_name}_error.txt', 'w') as notefile:
                    notefile.write(f'{dataset_name}\ntarget_name {target_name}\n')
                    notefile.write(f'entry {entry_index}\n')
                    notefile.write(f'error message:\n{msg}')
                # remove this folder
                os.chdir('..')
                shutil.rmtree(target_name)
                print("Error: " + msg.replace('\n', ';'))
            else:
                # write conf.pdb file
                fbmol = Molecule(f'input.mol2')
                fbmol.write(f'conf.pdb')
                # write vdata.txt for this target
                create_vdata_txt('conf.pdb', data)
                print(f'Success')
                # finish this target
                target_names.append(target_name)
                os.chdir('..')
            target_idx += 1
        except: 
            if not os.path.exists('../error_mol2s'):
                os.mkdir('../error_mol2s')
            err_msg = 'cmiles can not load the molecule. The molecule may have elements that are unacceptable. (e.g. Si)'            
            with open(f'../error_mol2s/{target_name}_error.txt', 'w') as notefile:
                notefile.write(f'{dataset_name}\ntarget_name {target_name}\n')
                notefile.write(f'entry {entry_index}\n')
                notefile.write(f'error message:\n{errmsg}')
            # remove this folder
            os.chdir('..')
            shutil.rmtree(target_name)
            print("Error: " + errmsg.replace('\n', ';'))

    # write targets.{dataset_name}.in file
    target_in_fnm = f"targets.{dataset_name.replace(' ', '_')}.in"
    with open(target_in_fnm, 'w') as outfile:
        for target_name in target_names:
            outfile.write(target_in_str.format(name=target_name))
    print(f"Successfull generated {len(target_names)} targets.")
    print(f"You can copy contents in {target_in_fnm} to your ForceBalance input file.")
    os.chdir('..')


def get_int_fmt_string(n):
    # count number of digits needed
    count = 0
    while n > 0:
        n //= 10
        count += 1
    return f"{{:0{count}d}}"

def check_connectivity(filename):
    """ Check if the connectivity in the molecule file is consistent with geometry
    Using force balance Molecule.build_bonds() for this first draft
    This can be improved by OpenEye or other methods
    """
    fbmol = Molecule(filename)
    orig_bonds = set(fbmol.bonds)
    # for b1, b2 in fbmol.bonds:
    #     bond = (b1, b2) if b1 < b2 else (b2, b1)
    #     orig_bonds.add(bond)
    fbmol.build_bonds()
    new_bonds = set(fbmol.bonds)
    # for b1, b2 in fbmol.bonds:
    #     bond = (b1, b2) if b1 < b2 else (b2, b1)
    #     new_bonds.add(bond)
    return orig_bonds == new_bonds

def test_ff_mol2(test_ff, mol2_fnm):
    """
    Test creating system with mol2 file
    """
    from openforcefield.topology import Molecule as Off_Molecule
    from openforcefield.topology import Topology as Off_Topology
    try:
        off_molecule = Off_Molecule.from_file(mol2_fnm)
        off_topology = Off_Topology.from_molecules(off_molecule)
        test_ff.create_openmm_system(off_topology)
        molecule_labels = test_ff.label_molecules(off_topology)[0]
    except Exception as e:
        return False, str(e), None
    return True, '', molecule_labels

def create_vdata_txt(pdb_fnm, data):
    # check gradient
    gradient = data['gradient']
    if np.abs(gradient).max() > 1e-3:
        print("Warning! Max gradient greater than 1e-3")
    # get the openmm list of masses for the molecule, to be consistent with ForceBalance
    pdb = app.PDBFile(pdb_fnm)
    mass_array = np.array([a.element.mass.value_in_unit(unit.dalton) for a in pdb.topology.atoms()])
    # compute mass-weighted hessian
    invert_sqrt_mass_array_repeat = 1.0 / np.sqrt(mass_array.repeat(3)) # repeat for x, y, z coords
    hessian = data['hessian']
    mass_weighted_hessian = hessian * invert_sqrt_mass_array_repeat[:, np.newaxis] * invert_sqrt_mass_array_repeat[np.newaxis, :]
    mass_weighted_hessian *= 937583.07 # convert units from Eh / bohr^2 * dalton to 10^24 s^-2
    # perform normal mode analysis
    freqs, normal_modes = compute_normal_modes(mass_weighted_hessian)
    # write vdata.txt
    with open('vdata.txt', 'w') as outfile:
        outfile.write(data['molecule'].to_string('xyz') + '\n')
        for freq, normol_mode in zip(freqs, normal_modes):
            outfile.write(f'{freq}\n')
            for nx, ny, nz in normol_mode:
                outfile.write(f'{nx:13.4f} {ny:13.4f} {nz:13.4f}\n')
            outfile.write('\n')

def compute_normal_modes(mass_weighted_hessian):
    # 1. diagonalize the hessian matrix
    eigvals, eigvecs = np.linalg.eigh(mass_weighted_hessian)
    # 2. convert eigenvalues to frequencies
    coef = 0.5 / np.pi * 33.3564095 # 10^12 Hz => cm-1
    negatives = (eigvals >= 0).astype(int) * 2 - 1 # record the negative ones
    freqs = np.sqrt(np.abs(eigvals)) * coef * negatives
    # 3. convert eigenvectors to normal modes
    noa = int(len(mass_weighted_hessian)/3)
    # scale the eigenvectors by mass^-1/2 (This is not needed)
    #normal_modes = (eigvecs.reshape(noa, -1) / np.sqrt(mass_array)[:, np.newaxis]).reshape(noa*3, noa*3)
    normal_modes = eigvecs
    # re-arange to row index and shape
    normal_modes = normal_modes.T.reshape(noa*3, noa, 3)
    # nomalize the scaled normal modes
    # norm_factors = np.sqrt(np.sum(np.square(normal_modes), axis=(1,2)))
    # normal_modes /= norm_factors[:, np.newaxis, np.newaxis]
    # step 5: remove the 6 freqs with smallest abs value and corresponding normal modes
    n_remove = 5 if noa == 2 else 6
    larger_freq_idxs = np.sort(np.argpartition(np.abs(freqs), n_remove)[n_remove:])
    freqs = freqs[larger_freq_idxs]
    normal_modes = normal_modes[larger_freq_idxs]
    return freqs, normal_modes


def main():
    import argparse
    import pickle
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset", help='Name of the OptimizationDataset on QCFractal')
    parser.add_argument("-l", "--load_pickle", help='Load downloaded torsiondrive data from pickle file')
    parser.add_argument("-t", "--test_ff_fnm", required=True, help="Provide an offxml for testing the molecules created, skip the ones that failed")
    args = parser.parse_args()

    # step 1: download dataset from server
    if args.load_pickle:
        with open(args.load_pickle, 'rb') as pfile:
            hessian_data = pickle.load(pfile)
            print(f"{len(hessian_data)} hessian data entries loaded from {args.load_pickle}")
    else:
        hessian_data = download_hessian_data(args.dataset)
    # step 2: filter the dataset, keep one lowest energy for each molecule
    hessian_data = filter_hessian_data(hessian_data)

    # create a ForceField object for testing
    # test_ff = ForceField(args.test_ff_fnm)
    test_ff = ForceField(args.test_ff_fnm, allow_cosmetic_attributes=True)
    # step 3: generate one target for each data entry
    make_vib_freq_target(args.dataset, hessian_data, test_ff=test_ff)

if __name__ == '__main__':
    main()
