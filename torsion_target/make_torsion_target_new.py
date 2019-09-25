#!/usr/bin/env python

import os
import shutil
import pickle
import numpy as np
import copy
import json

import cmiles
from openeye import oechem
import qcportal as ptl
from forcebalance.molecule import Molecule
from openforcefield.typing.engines.smirnoff import ForceField

client = ptl.FractalClient('https://api.qcarchive.molssi.org:443/')
ofs = oechem.oemolostream()


def download_torsiondrive_data(dataset_name):
    """
    Download data from public server

    Parameters
    ----------
    dataset_name: str
        example: "SMIRNOFF Coverage Torsion Set 1"

    Returns
    -------
    torsiondrive_data: dict
        {
            '[H:4][CH:3]1C(=O)N[C@@:2]1(C)[c:1]2nccs2': {
                'initial_molecules': [
                    <Molecule(name='C7H8N2OS' formula='C7H8N2OS' hash='3a1da45')>,
                    <Molecule(name='C7H8N2OS' formula='C7H8N2OS' hash='30b8ded')>,
                ],
                'final_molecules': {
                    (-135,): <Molecule(name='C7H8N2OS' formula='C7H8N2OS' hash='dcff9dd')>,
                    (-150,): <Molecule(name='C7H8N2OS' formula='C7H8N2OS' hash='c166072')>,
                    (-120,): <Molecule(name='C7H8N2OS' formula='C7H8N2OS' hash='f7b637e')>,
                    ...
                },
                'final_energies': {
                    (-135,): -854.5232944059941,
                    (-150,): -854.5209682017379,
                    (-120,): -854.523412459059,
                    ...
                },
                'final_gradients': {
                    (-135,): np.array([2.27e-5, -4.69e-5, ...]),
                    (-150,): np.array([0.33e-5, -5.93e-5, ...]),
                    (-120,): np.array([4.21e-5, -1.22e-5, ...]),
                },
                'keywords': {
                    'dihedrals': [(2, 5, 4, 14)],
                    'grid_spacing': [15],
                    'dihedral_ranges': None,
                    'energy_decrease_thresh': None,
                    'energy_upper_limit': 0.05,
                },
                'attributes': {
                    'canonical_explicit_hydrogen_smiles': '[H]C([H])([H])OC([H])(O[H])OC([H])([H])[H]',
                    'canonical_isomeric_explicit_hydrogen_mapped_smiles': '[H:7][C:1]([H:8])([H:9])[O:5][C:3]([H:13])([O:4][H:14])[O:6][C:2]([H:10])([H:11])[H:12]',
                    'canonical_isomeric_explicit_hydrogen_smiles': '[H]C([H])([H])OC([H])(O[H])OC([H])([H])[H]',
                    'canonical_isomeric_smiles': 'COC(O)OC',
                    'canonical_smiles': 'COC(O)OC',
                    'inchi_key': 'IIGJYLXJNYBXEO-UHFFFAOYSA-N',
                    'molecular_formula': 'C3H8O3',
                    'provenance': 'cmiles_v0.1.5_openeye_2019.Apr.2',
                    'standard_inchi': 'InChI=1S/C3H8O3/c1-5-3(4)6-2/h3-4H,1-2H3',
                    'unique_protomer_representation': 'COC(O)OC',
                    'unique_tautomer_representation': 'COC(O)OC',
                },
            },
            ...
        }
    """
    # load dataset from public qcfractal server
    ds = client.get_collection("TorsionDriveDataset", dataset_name)
    spec_name = ds.list_specifications().index[0]
    print(f"Loading TorsionDrive Scans from [ {dataset_name} ] spec [{spec_name}]")
    print(f"Found {len(ds.df)} data entries")
    # load torsiondrive record ids from the dataset
    map_record_id_entry_index = {}
    for entry_index in ds.df.index:
        data_entry = ds.get_entry(entry_index)
        td_record_id = data_entry.object_map[spec_name]
        map_record_id_entry_index[td_record_id] = entry_index, data_entry.attributes
    print(f"Found {len(map_record_id_entry_index)} torsiondrive records")
    # query all torsiondrive records at the same time
    td_record_ids = list(map_record_id_entry_index.keys())
    torsiondrive_data = {}
    for i, td_record in enumerate(client.query_procedures(id=td_record_ids), 1):
        entry_index, attributes = map_record_id_entry_index[td_record.id]
        print(f"{i:5d} : {entry_index:50s} status {td_record.status}")
        if td_record.status == 'COMPLETE':
            torsiondrive_data[entry_index] = {
                'initial_molecules': client.query_molecules(td_record.initial_molecule),
                'final_molecules': td_record.get_final_molecules(),
                'final_energies': td_record.get_final_energies(),
                'final_gradients': {gid: np.array(res.return_result) for gid, res in td_record.get_final_results().items()},
                'keywords': td_record.keywords.dict(),
                'attributes': attributes,
            }
    print(f'Downloaded torsion drive data for {len(torsiondrive_data)} completed entries')
    # save as pickle file
    with open('torsiondrive_data.pickle', 'wb') as pfile:
        pickle.dump(torsiondrive_data, pfile)
    return torsiondrive_data

target_in_str = '''
$target
name {name}
type TorsionProfile_SMIRNOFF
mol2 input.mol2
pdb conf.pdb
coords scan.xyz
writelevel 2
attenuate
energy_denom 1.0
energy_upper 5.0
remote 1
$end
'''

def make_torsiondrive_target(dataset_name, torsiondrive_data, test_ff=None):
    """
    Make a folder of ForceBalance targets from the torsiondrive data
    """
    target_name_prefix = 'td_' + dataset_name.replace(' ', '_')
    # create new targets folder
    if os.path.exists('targets'):
        shutil.rmtree('targets')
    os.mkdir('targets')
    os.chdir('targets')
    # write each entry as an individual target
    target_idx = 0
    n_targets = len(torsiondrive_data)
    idx_fmt_string = get_int_fmt_string(n_targets)
    target_names = []
    print(f"Generating {n_targets} targets")
    for entry_index, td_data in torsiondrive_data.items():
        # pick a single initial molecule
        qcmol = td_data['initial_molecules'][0]
        # get mol_formula
        mol_formula = qcmol.get_molecular_formula()
        # create target folder
        target_idx_str = idx_fmt_string.format(target_idx)
        target_name = f"{target_name_prefix}_{target_idx_str}_{mol_formula}"
        print(f"{target_idx}: {target_name}")
        os.mkdir(target_name)
        os.chdir(target_name)
        # write a note
        with open('note.txt', 'w') as notefile:
            notefile.write(f'Target generated from dataset {dataset_name}, entry {entry_index}')
        # write input.mol2 file
        qcjson_mol = qcmol.dict(encoding='json')
        oemol = cmiles.utils.load_molecule(qcjson_mol)
        ofs.open(f'input.mol2')
        oechem.OEWriteMolecule(ofs, oemol)
        ofs.close()
        # test mol2 file
        success = True
        if test_ff != None:
            success, msg, molecule_labels = test_ff_mol2(test_ff, 'input.mol2')
        # check if the torsion scan contains one or more conformers forming strong internal H bonds
        if success:
            # write conf.pdb file
            fbmol = Molecule(f'input.mol2')
            fbmol.write(f'conf.pdb')
            # list of grid ids sorted
            sorted_grid_ids = sorted(td_data['final_molecules'].keys())
            # write scan.xyz and qdata.txt files
            target_mol = Molecule()
            target_mol.elem = fbmol.elem
            target_mol.xyzs = []
            target_mol.qm_energies = []
            target_mol.qm_grads = []
            for grid_id in sorted_grid_ids:
                grid_qc_mol = td_data['final_molecules'][grid_id]
                # convert geometry unit Bohr -> Angstrom
                geo = grid_qc_mol.geometry * 0.529177
                target_mol.xyzs.append(geo)
                # add energy and gradient
                target_mol.qm_energies.append(td_data['final_energies'][grid_id])
                target_mol.qm_grads.append(td_data['final_gradients'][grid_id])
            target_mol.write('scan.xyz')
            target_mol.write('qdata.txt')

            no_hbonds = check_Hbond(scan_fnm ='scan.xyz', top_fnm='input.mol2')
            if not no_hbonds:
                success = False
                msg = 'One or more internal H bonds exist.'
        if not success:
            if not os.path.exists('../error_mol2s'):
                os.mkdir('../error_mol2s')
            shutil.move(f'input.mol2', f'../error_mol2s/{target_name}.mol2')
            with open(f'../error_mol2s/{target_name}_error.txt', 'w') as notefile:
                notefile.write(f'{dataset_name}\ntarget_name {target_name}\n')
                notefile.write(f'entry {entry_index}\ntd_keywords {td_data["keywords"]}\n')
                notefile.write(f'error message:\n{msg}')
            # remove this folder
            os.chdir('..')
            shutil.rmtree(target_name)
        else:
            # pick metadata to write into the metadata.json file
            metadata = copy.deepcopy(td_data['keywords'])
            metadata['dataset_name'] = dataset_name
            metadata['entry_label'] = entry_index
            metadata['canonical_smiles'] = td_data['attributes'].get('canonical_smiles', 'unknown')
            metadata['attributes'] = td_data['attributes']
            metadata['torsion_grid_ids'] = sorted_grid_ids
            # find SMIRKs for torsion being scaned if test_ff is provided
            if test_ff:
                metadata['smirks'] = []
                metadata['smirks_ids'] = []
                for torsion_indices in td_data['keywords']['dihedrals']:
                    param = molecule_labels['ProperTorsions'][tuple(torsion_indices)]
                    metadata['smirks'].append(param.smirks)
                    metadata['smirks_ids'].append(param.id)
            with open('metadata.json', 'w') as jsonfile:
                json.dump(metadata, jsonfile, indent=2)
            # finish this target
            target_names.append(target_name)
            os.chdir('..')
        target_idx += 1

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

def check_Hbond(scan_fnm, top_fnm=None):
    """
    Check if the torsion scan contains conformers with internal hydrogen bonds
    """
    import mdtraj as md
    traj = md.load(scan_fnm, top=top_fnm)
    hbonds = md.baker_hubbard(traj)
    if len(hbonds) == 0:
        return True
    else:
        return False

def main():
    import argparse
    import pickle
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset", help='Name of the OptimizationDataset on QCFractal')
    parser.add_argument("-l", "--load_pickle", help='Load downloaded torsiondrive data from pickle file')
    parser.add_argument("-t", "--test_ff_fnm", required=True, help="Provide an offxml for testing the molecules created, skip the ones that failed")
    args = parser.parse_args()

    if args.load_pickle:
        with open(args.load_pickle, 'rb') as pfile:
            torsiondrive_data = pickle.load(pfile)
    else:
        torsiondrive_data = download_torsiondrive_data(args.dataset)

    # require the test_ff for metadata
    test_ff = ForceField(args.test_ff_fnm)
    make_torsiondrive_target(args.dataset, torsiondrive_data, test_ff=test_ff)

if __name__ == '__main__':
    main()
