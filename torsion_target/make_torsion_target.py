#!/usr/bin/env python

import os
import shutil
import pickle
import numpy as np

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
                }
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
        map_record_id_entry_index[td_record_id] = entry_index
    print(f"Found {len(map_record_id_entry_index)} torsiondrive records")
    # query all torsiondrive records at the same time
    td_record_ids = list(map_record_id_entry_index.keys())
    torsiondrive_data = {}
    for i, td_record in enumerate(client.query_procedures(id=td_record_ids), 1):
        entry_index = map_record_id_entry_index[td_record.id]
        print(f"{i:5d} : {entry_index:30s} status {td_record.status}")
        if td_record.status == 'COMPLETE':
            torsiondrive_data[entry_index] = {
                'initial_molecules': client.query_molecules(td_record.initial_molecule),
                'final_molecules': td_record.get_final_molecules(),
                'final_energies': td_record.get_final_energies(),
                'final_gradients': {gid: np.array(res.return_result) for gid, res in td_record.get_final_results().items()}
            }
    print(f'Downloaded torsion drive data for {len(torsiondrive_data)} completed entries')
    # save as pickle file
    with open('torsiondrive_data.pickle', 'wb') as pfile:
        pickle.dump(torsiondrive_data, pfile)
    return torsiondrive_data

target_in_str = '''
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
        qcjson_mol = qcmol.json_dict()
        oemol = cmiles.utils.load_molecule(qcjson_mol)
        ofs.open(f'input.mol2')
        oechem.OEWriteMolecule(ofs, oemol)
        ofs.close()
        # test mol2 file
        success = True
        if test_ff != None:
            success, msg = test_ff_mol2(test_ff, 'input.mol2')
        if not success:
            if not os.path.exists('../error_mol2s'):
                os.mkdir('../error_mol2s')
            shutil.move(f'input.mol2', f'../error_mol2s/{target_name}.mol2')
            with open(f'../error_mol2s/{target_name}_error.txt', 'w') as notefile:
                notefile.write(f'{dataset_name}\nentry {entry_index}\ntarget_name {target_name}\n')
                notefile.write(msg)
            # remove this folder
            os.chdir('..')
            shutil.rmtree(target_name)
        else:
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
            target_names.append(target_name)
            # finish this target
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
    except Exception as e:
        return False, str(e)
    return True, ''


def main():
    import argparse
    import pickle
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset", help='Name of the OptimizationDataset on QCFractal')
    parser.add_argument("-l", "--load_pickle", help='Load downloaded torsiondrive data from pickle file')
    parser.add_argument("-t", "--test_ff_fnm", help="Provide an offxml for testing the molecules created, skip the ones that failed")
    args = parser.parse_args()

    if args.load_pickle:
        with open(args.load_pickle, 'rb') as pfile:
            torsiondrive_data = pickle.load(pfile)
    else:
        torsiondrive_data = download_torsiondrive_data(args.dataset)

    test_ff = ForceField(args.test_ff_fnm) if args.test_ff_fnm else None
    make_torsiondrive_target(args.dataset, torsiondrive_data, test_ff=test_ff)

if __name__ == '__main__':
    main()
