#!/usr/bin/env python

import os
import cmiles
from openeye import oechem
import qcportal as ptl
client = ptl.FractalClient('https://api.qcarchive.molssi.org:443/')
from openforcefield.topology import Molecule
from forcebalance.molecule import Molecule

ofs = oechem.oemolostream()

global_opts = """
$global
bond_denom 0.01
angle_denom 0.03
dihedral_denom 0.2
improper_denom 0.2
$end
"""

target_template = """
$system
name {name}
geometry {name}.xyz
topology {name}.pdb
mol2 {name}.mol2
$end
"""

target_opts = """
$target
name {name}
type OptGeoTarget_SMIRNOFF
weight 1.0
writelevel 1
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
    for optRecord in client.query_procedures(id=optRecord_ids):
        m_index = map_optRecordId_mIndex[optRecord.id]
        if optRecord.final_molecule:
            final_molecules_ids_map[optRecord.final_molecule] = m_index
    # query all final molecules
    final_molecules = {}
    final_molecule_ids = list(final_molecules_ids_map.keys())
    for molecule in client.query_molecules(final_molecule_ids):
        m_index = final_molecules_ids_map[molecule.id]
        final_molecules[m_index] = molecule
    print(f"Loaded {len(final_molecules)} final molecules")
    return final_molecules

def get_int_fmt_string(n):
    # count number of digits needed
    count = 0
    while n > 0:
        n //= 10
        count += 1
    return f"{{:0{count}d}}"

def write_molecule_files(molecule, name):
    qcjson_mol = molecule.json_dict()
    oemol = cmiles.utils.load_molecule(qcjson_mol)
    # write the mol2 file and xyz file
    for ext in ['xyz', 'mol2']:
        ofs.open(f'{name}.{ext}')
        oechem.OEWriteMolecule(ofs, oemol)
        ofs.close()
    # use ForceBalance Molecule to generate pdb file
    # because the oechem.OEWriteMolecule will mess up with atom indices
    fbmol = Molecule(f'{name}.mol2')
    fbmol.write(f'{name}.pdb')


def make_optgeo_target(dataset_name):
    final_molecules = load_final_molecules(dataset_name)
    # get into targets folder
    if not os.path.exists('targets'):
        os.mkdir('targets')
    os.chdir('targets')
    # prepare output optget target folder
    target_name = 'optgeo_' + dataset_name.replace(' ', '_')
    target_folder = os.path.join(target_name)
    if not os.path.exists(target_folder):
        os.mkdir(target_folder)
    os.chdir(target_folder)
    print(f"Generating optgeo target in {target_folder}")
    # write notes
    fnotes = open('notes.txt', 'w')
    fnotes.write("Prepared by make_fb_optgeo_target.py\n")
    fnotes.write(f"Data loaded from {dataset_name}\n")
    with open("optgeo_options.txt", "w") as optfile:
        optfile.write(global_opts)
        target_index = 0
        # index format with correct number of leading 0s
        idx_fmt_string = get_int_fmt_string(len(final_molecules))
        for m_index, molecule in final_molecules.items():
            idx_str = idx_fmt_string.format(target_index)
            name = f"{idx_str}_{molecule.name}"
            # if molecule.molecular_charge != 0:
            #     fnotes.write(f'{name} : {m_index} skipped - charge {molecule.molecular_charge}\n')
            #     continue
            # write molecule files
            write_molecule_files(molecule, name)
            # write $system block in optgeo_options.txt
            optfile.write(target_template.format(name=name))
            # write the original m_index
            fnotes.write(f'{name} : {m_index}\n')
            target_index += 1
    fnotes.close()
    os.chdir('..')
    # write a targets.in file
    with open('targets.in', 'w') as outfile:
        outfile.write(target_opts.format(name=target_name))
    os.chdir('..')
    print(f"Targets generation finished!")
    print(f"You can copy contents in {os.path.join('targets', 'targets.in')} to your ForceBalance input file.")



def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset", help='Name of the OptimizationDataset on QCFractal')
    args = parser.parse_args()

    make_optgeo_target(args.dataset)

if __name__ == '__main__':
    main()
