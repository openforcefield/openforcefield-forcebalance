#!/usr/bin/env python


# This script downloads data from an OptimizationDataset on QCArchieve server and make a list of ForceBalance AbInitio targets
# Step 1: Download data, including final molecule (geometry), energy
# Step 2: Combine conformers into molecules (each entry correspond to one conformer)
# Step 3: Create one Abinitio target for each molecule with more than 2 conformers. The QM vs MM relative energies will be used in the evaluation.

# This script is currently not used in fitting, but mainly for the purpose of benchmarking

import os
import shutil
import pickle
import json
import copy

import cmiles
from openeye import oechem
import qcportal as ptl
from forcebalance.molecule import Molecule, Elements, bohr2ang

client = ptl.FractalClient('https://api.qcarchive.molssi.org:443/')

ofs = oechem.oemolostream()

# we used energy_rms_override here mainly for benchmarking
# otherwise some conformers have very similar energy will cause a extremely small rms thus huge objective
target_opts = """
$target
name {name}
type AbInitio_SMIRNOFF
mol2 input.mol2
pdb conf.pdb
coords coords.xyz
writelevel 2
energy_rms_override 10.0
force 0
openmm_platform Reference
remote 1
$end
"""

def load_final_molecule_data(dataset_name):
    """ Download final molecule and energy from OptimizationDataset """
    # load dataset from public qcfractal server
    ds = client.get_collection("OptimizationDataset", dataset_name)
    spec_name = ds.list_specifications().index[0]
    print(f"Loading final molecules from [ {dataset_name} ] spec [{spec_name}]")
    print(f"Found {len(ds.df)} data entries")
    # load optimization record ids from the dataset
    map_optRecordId_mIndex = {}
    for m_index in ds.df.index:
        data_entry = ds.get_entry(m_index)
        optRecordId = data_entry.object_map[spec_name]
        map_optRecordId_mIndex[optRecordId] = m_index
    print(f"Found {len(map_optRecordId_mIndex)} optimization records")
    # query all opt records at the same time
    optRecord_ids = list(map_optRecordId_mIndex.keys())
    final_molecules_ids_map = {}
    final_molecule_energies = {}
    # query procedures by chunks
    chunk_size = 1000
    n_chunks = int((len(optRecord_ids)+chunk_size-1) / chunk_size)
    for i_chunk in range(n_chunks):
        chunk_record_ids = optRecord_ids[i_chunk * chunk_size: (i_chunk+1)*chunk_size]
        for optRecord in client.query_procedures(id=chunk_record_ids):
            m_index = map_optRecordId_mIndex[optRecord.id]
            if optRecord.final_molecule:
                final_molecules_ids_map[optRecord.final_molecule] = m_index
                final_molecule_energies[optRecord.final_molecule] = optRecord.energies[-1]
    # query all final molecules
    final_molecule_data = {}
    final_molecule_ids = list(final_molecules_ids_map.keys())
    print(f"Found {len(final_molecule_ids)} final_molecule_ids")
    # load by chunks of 1000 (limit on server)
    n_chunks = int((len(final_molecule_ids)+chunk_size-1) / chunk_size)
    for i_chunk in range(n_chunks):
        chunk_molecule_ids = final_molecule_ids[i_chunk * chunk_size: (i_chunk+1)*chunk_size]
        for molecule in client.query_molecules(chunk_molecule_ids):
            m_index = final_molecules_ids_map[molecule.id]
            energy = final_molecule_energies[molecule.id]
            final_molecule_data[m_index] = (molecule, energy)
    print(f"Loaded {len(final_molecule_data)} final molecules")
    # save as pickle file
    with open('final_molecule_data.pickle', 'wb') as pfile:
        pickle.dump(final_molecule_data, pfile)
    return final_molecule_data

def get_int_fmt_string(n):
    # count number of digits needed
    count = 0
    while n > 0:
        n //= 10
        count += 1
    return f"{{:0{count}d}}"

def group_label_into_molecules(conformer_entry_labels):
    """
    Take a list of conformer entry labels, use the pattern XXX-0, XXX-1 to determine which molecule it belongs to
    Return a dictionary of {molecule_label: [conformer_labels]}
    """
    molecule_conformer_dict = {}
    for l in conformer_entry_labels:
        molecule_label = l.rsplit('-', maxsplit=1)[0]
        molecule_conformer_dict.setdefault(molecule_label, [])
        molecule_conformer_dict[molecule_label].append(l)
    return molecule_conformer_dict

def make_abinitio_targets(dataset_name, final_molecule_data, test_ff_fnm=None):
    # create the ff for testing
    if test_ff_fnm != None:
        from openforcefield.typing.engines.smirnoff import ForceField
        test_ff = ForceField(test_ff_fnm, allow_cosmetic_attributes=True)
    else:
        test_ff = None
    # group conformer labels into molecules
    molecule_conformer_dict = group_label_into_molecules(sorted(final_molecule_data.keys()))
    n_molecules = len(molecule_conformer_dict)
    if n_molecules == 0:
        print("No molecules found")
        return
    # create and cd into targets folder
    if not os.path.exists('targets'):
        os.mkdir('targets')
    os.chdir('targets')
    idx_fmt_string = get_int_fmt_string(n_molecules)
    i_target = 0
    n_success = 0
    target_names = []
    for molecule_label, conformer_entry_labels in molecule_conformer_dict.items():
        if len(conformer_entry_labels) < 3:
            print(f"molecule {molecule_label} skipped. Reason: Length of conformer_entry_labels {len(conformer_entry_labels)} < 3")
            continue
        #  create one target for each molecule
        molecule_conformer_data = {l: final_molecule_data[l] for l in conformer_entry_labels}
        target_name = 'abinitio_' + dataset_name.replace(' ', '_') + '-' + idx_fmt_string.format(i_target)
        success = create_one_abinitio_target(target_name, molecule_conformer_data)
        if success:
            n_success += 1
            target_names.append(target_name)
        else:
            # move failed targets into a separate folder
            if not os.path.exists('error_targets'):
                os.mkdir('error_targets')
            if os.path.isdir(f'{target_name}'):
                shutil.move(f'{target_name}', f'error_targets/{target_name}')
            print(f"Create target {target_name} failed, check targets/error_targets/{target_name} for details")
        i_target += 1
    print(f"Successfully created targets for {n_success} molecules")
    # write a targets.in file
    target_in_fnm = f"targets.in.abinitio.{dataset_name.replace(' ', '_')}"
    with open(target_in_fnm, 'w') as outfile:
        for target_name in target_names:
            outfile.write(target_opts.format(name=target_name))
    os.chdir('..')
    print(f"Targets generation finished!")
    print(f"You can copy contents in {os.path.join('targets', target_in_fnm)} to your ForceBalance input file.")


def create_one_abinitio_target(target_name, molecule_conformer_data, test_ff=None):
    """ generate a single target folder with data provided
    moledules_data: [str: ] = {molecule_index : Molecule}
    """
    assert len(molecule_conformer_data) > 2
    # prepare target folder
    target_folder = target_name
    if not os.path.exists(target_folder):
        os.mkdir(target_folder)
    os.chdir(target_folder)
    print(f"Generating abinitio target in {target_folder}")
    # prepare notes
    fnotes = open('notes.txt', 'w')
    fnotes.write("Prepared by make_fb_abinitio_target.py\n")
    fnotes.write(f"Relative energies of {len(molecule_conformer_data)} conformers\n")
    # filter molecule list by a few checks
    molecule_data_list = []
    success = True
    for conformer_label, mdata in molecule_conformer_data.items():
        m, e = mdata
        success, err_msg = check_molecule(m, test_ff)
        if success:
            molecule_data_list.append(mdata)
            fnotes.write(f'{conformer_label}: molecule_id {m.id} | SUCCESS\n')
        else:
            fnotes.write(f'{conformer_label}: molecule_id {m.id} | ERROR: {err_msg}\n')
    if len(molecule_data_list) < 2:
        fnotes.write("Remaining conformers < 2 after filerting, skipping this target\n")
        success = False
    if success:
        # write files for remaining molecules
        write_molecule_files(molecule_data_list)
    # finish
    fnotes.close()
    os.chdir("..")
    return success

def check_molecule(molecule, test_ff=None):
    """ run a few checks for a QCElemental Molecule """
    import tempfile
    qcjson_mol = molecule.dict(encoding='json')
    oemol = cmiles.utils.load_molecule(qcjson_mol)
    success = True
    err_msg = ""
    cwd = os.getcwd()
    # write a test.mol2 file in a temp dir for checking
    with tempfile.TemporaryDirectory() as tmpdirname:
        os.chdir(tmpdirname)
        ofs.open('test.mol2')
        oechem.OEWriteMolecule(ofs, oemol)
        ofs.close()
        # test if bonds changed
        if not check_connectivity('test.mol2'):
            success = False
            err_msg = "Bonds changed after rebuild"
        # test if can be created by the test_ff
        if success == True and test_ff != None:
            from openforcefield.topology import Molecule as Off_Molecule
            from openforcefield.topology import Topology as Off_Topology
            try:
                off_molecule = Off_Molecule.from_file(f'{name}.mol2')
                off_topology = Off_Topology.from_molecules(off_molecule)
                test_ff.create_openmm_system(off_topology)
            except Exception as e:
                success = False
                err_msg = str(e)
        # test if this molecule has hydrogen bonds
        if not check_hbond('test.mol2'):
            success = False
            err_msg = 'One or more hydrogen bond found'
    # go back to orig dir
    os.chdir(cwd)
    return success, err_msg

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

def check_hbond(mol2_fnm):
    import mdtraj as md
    traj = md.load(mol2_fnm)
    hbonds = md.baker_hubbard(traj)
    return len(hbonds) == 0

def write_molecule_files(molecule_data_list):
    molecule, e0 = molecule_data_list[0]
    qcjson_mol = molecule.dict(encoding='json')
    oemol = cmiles.utils.load_molecule(qcjson_mol)
    # write the mol2 file using oechem
    ofs.open(f'input.mol2')
    oechem.OEWriteMolecule(ofs, oemol)
    ofs.close()
    # write the pdb file using ForceBalance Molecule
    fbmol = Molecule(f'input.mol2')
    fbmol.write(f'conf.pdb')
    # write xyz file using a new ForceBalance Molecule object
    m = Molecule()
    m.elem = [Elements[i] for i in molecule.atomic_numbers]
    m.xyzs = []
    m.qm_energies = []
    for mol, e in molecule_data_list:
        m.xyzs.append(mol.geometry * bohr2ang)
        m.qm_energies.append(e)
    m.write("coords.xyz")
    # write qdata.txt file with coords and energies
    m.write('qdata.txt')

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset", help='Name of the OptimizationDataset on QCFractal')
    parser.add_argument("-l", "--load_pickle", help='Load downloaded torsiondrive data from pickle file')
    parser.add_argument("-s", "--size", default=50, type=int, help="Size of each target, equal to number of optimized geometries in each target")
    parser.add_argument("-t", "--test_ff_fnm", help="Provide an offxml for testing the molecules created, skip the ones that failed")
    args = parser.parse_args()

    if args.load_pickle:
        with open(args.load_pickle, 'rb') as pfile:
            final_molecule_data = pickle.load(pfile)
    else:
        final_molecule_data = load_final_molecule_data(args.dataset)

    make_abinitio_targets(args.dataset, final_molecule_data, test_ff_fnm=args.test_ff_fnm)

if __name__ == '__main__':
    main()
