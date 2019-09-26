#!/usr/bin/env python

import os
import shutil
import pickle

import cmiles
from openeye import oechem
import qcportal as ptl
client = ptl.FractalClient('https://api.qcarchive.molssi.org:443/')
from forcebalance.molecule import Molecule

ofs = oechem.oemolostream()

# NOTE by QYD:
# 1. dihedral_denom is 0 for now because equilibrium dihedral angle
#    does not have a well-defined gradient against Torsion k when k is near 0
#    i.e. k=+0.001  =>  equilibrium dihedral angle = 0.0
#         k=-0.001  =>  equilibrium dihedral angle = 180.0
# 2. The bond_denom and angle_denom is estimated to have similar energy impact
#    based on bond_k ~ 1000/ang^2  angle_k ~100/radian^2
global_opts = """
$global
bond_denom 0.05
angle_denom 8
dihedral_denom 0
improper_denom 20
$end
"""

system_template = """
$system
name {name}
geometry {name}.xyz
topology {name}.pdb
mol2 {name}.mol2
$end
"""

# Note: weight is adjusted lower to prevent dominating the total objective
target_opts = """
$target
name {name}
type OptGeoTarget_SMIRNOFF
weight 0.1
writelevel 1
remote 1
$end
"""

def load_final_molecules(dataset_name):
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
    # query procedures by chunks
    chunk_size = 1000
    n_chunks = int((len(optRecord_ids)+chunk_size-1) / chunk_size)
    for i_chunk in range(n_chunks):
        chunk_record_ids = optRecord_ids[i_chunk * chunk_size: (i_chunk+1)*chunk_size]
        for optRecord in client.query_procedures(id=chunk_record_ids):
            m_index = map_optRecordId_mIndex[optRecord.id]
            if optRecord.final_molecule:
                final_molecules_ids_map[optRecord.final_molecule] = m_index
    # query all final molecules
    final_molecules = {}
    final_molecule_ids = list(final_molecules_ids_map.keys())
    print(f"Found {len(final_molecule_ids)} final_molecule_ids")
    # load by chunks of 1000 (limit on server)
    n_chunks = int((len(final_molecule_ids)+chunk_size-1) / chunk_size)
    for i_chunk in range(n_chunks):
        chunk_molecule_ids = final_molecule_ids[i_chunk * chunk_size: (i_chunk+1)*chunk_size]
        for molecule in client.query_molecules(chunk_molecule_ids):
            m_index = final_molecules_ids_map[molecule.id]
            final_molecules[m_index] = molecule
    print(f"Loaded {len(final_molecules)} final molecules")
    # save as pickle file
    with open('final_molecules.pickle', 'wb') as pfile:
        pickle.dump(final_molecules, pfile)
    return final_molecules

def get_int_fmt_string(n):
    # count number of digits needed
    count = 0
    while n > 0:
        n //= 10
        count += 1
    return f"{{:0{count}d}}"

def write_molecule_files(molecule, name, test_ff=None):
    qcjson_mol = molecule.dict(encoding='json')
    oemol = cmiles.utils.load_molecule(qcjson_mol)
    # write the mol2 file and xyz file
    for ext in ['xyz', 'mol2']:
        ofs.open(f'{name}.{ext}')
        oechem.OEWriteMolecule(ofs, oemol)
        ofs.close()
    success = True
    err_msg = ""
    # test if bonds changed
    if not check_connectivity(f'{name}.mol2'):
        success = False
        err_msg = "Bonds changed after rebuild"
        if not os.path.exists('../err_bonds_changed'):
            os.mkdir('../err_bonds_changed')
        shutil.move(f'{name}.mol2', f'../err_bonds_changed/{name}.mol2')
        os.remove(f'{name}.xyz')
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
            if not os.path.exists('../error_mol2s'):
                os.mkdir('../error_mol2s')
            shutil.move(f'{name}.mol2', f'../error_mol2s/{name}.mol2')
            os.remove(f'{name}.xyz')
    if success == True:
        # use ForceBalance Molecule to generate pdb file
        # because the oechem.OEWriteMolecule will mess up with atom indices
        fbmol = Molecule(f'{name}.mol2')
        fbmol.write(f'{name}.pdb')
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

def make_optgeo_target(dataset_name, final_molecules, size=None, test_ff_fnm=None):
    # create the ff for testing
    if test_ff_fnm != None:
        from openforcefield.typing.engines.smirnoff import ForceField
        test_ff = ForceField(test_ff_fnm, allow_cosmetic_attributes=True)
    else:
        test_ff = None
    n_molecules = len(final_molecules)
    if n_molecules == 0:
        print("No molecules found")
        return
    # create and cd into targets folder
    if not os.path.exists('targets'):
        os.mkdir('targets')
    os.chdir('targets')
    if size is None or n_molecules <= size:
        # put all molecules into one target
        target_name = 'optgeo_' + dataset_name.replace(' ', '_')
        n_success = create_target(target_name, final_molecules)
        target_names = [target_name]
    else:
        assert size > 0 and isinstance(size, int), 'size should be positive int'
        n_groups = (n_molecules + size - 1) // size
        # put molecules into separate targets
        all_molelcule_idxs = list(final_molecules.keys())
        target_names = []
        idx_fmt_string = get_int_fmt_string(n_groups)
        n_success = 0
        for i_g in range(n_groups):
            group_molecule_idxs = all_molelcule_idxs[i_g*size : (i_g+1)*size]
            molecules_data = {m_index: final_molecules[m_index] for m_index in group_molecule_idxs}
            target_name = 'optgeo_' + dataset_name.replace(' ', '_') + '-' + idx_fmt_string.format(i_g)
            this_n_success = create_target(target_name, molecules_data, test_ff=test_ff, start_idx=i_g*size)
            if this_n_success > 0:
                target_names.append(target_name)
                n_success += this_n_success
    print(f"Successfully created targets with total {n_success} molecules")
    # write a targets.in file
    target_in_fnm = f"targets.in.{dataset_name.replace(' ', '_')}"
    with open(target_in_fnm, 'w') as outfile:
        for target_name in target_names:
            outfile.write(target_opts.format(name=target_name))
    os.chdir('..')
    print(f"Targets generation finished!")
    print(f"You can copy contents in {os.path.join('targets', target_in_fnm)} to your ForceBalance input file.")


def create_target(target_name, moledules_data, test_ff=None, start_idx=0):
    """ generate a single target folder with data provided
    moledules_data: [str: ] = {molecule_index : Molecule}
    """
    if len(moledules_data) == 0: return
    # prepare output optget target folder
    target_folder = os.path.join(target_name)
    if not os.path.exists(target_folder):
        os.mkdir(target_folder)
    os.chdir(target_folder)
    print(f"Generating optgeo target in {target_folder}")
    # write notes
    fnotes = open('notes.txt', 'w')
    fnotes.write("Prepared by make_fb_optgeo_target.py\n")
    n_success = 0
    with open("optgeo_options.txt", "w") as optfile:
        optfile.write(global_opts)
        target_index = start_idx
        # index format with correct number of leading 0s
        idx_fmt_string = get_int_fmt_string(len(moledules_data))
        for m_index, molecule in moledules_data.items():
            idx_str = idx_fmt_string.format(target_index)
            name = f"{idx_str}_{molecule.name}"
            # write molecule files
            success, err_msg = write_molecule_files(molecule, name, test_ff=test_ff)
            if not success:
                fnotes.write(f'{name} : {m_index} molecule_id {molecule.id} | ERROR: {err_msg}\n')
            else:
                # write $system block in optgeo_options.txt
                optfile.write(system_template.format(name=name))
                # write the original m_index
                fnotes.write(f'{name} : {m_index} molecule_id {molecule.id} | SUCCESS\n')
                n_success += 1
            target_index += 1
    fnotes.close()
    os.chdir("..")
    return n_success


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
            final_molecules = pickle.load(pfile)
    else:
        final_molecules = load_final_molecules(args.dataset)

    make_optgeo_target(args.dataset, final_molecules, size=args.size, test_ff_fnm=args.test_ff_fnm)

if __name__ == '__main__':
    main()
