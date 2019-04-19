#!/usr/bin/env python
# coding: utf-8
import os
import shutil
from openeye import oechem
ifs = oechem.oemolistream()
ofs = oechem.oemolostream()

# open input file that contains many molecules
ifs.open('OpenFF_references.sdf')

# get all molecules from sdf file
mol_list = []
for mol in ifs.GetOEGraphMols():
    # explicit declaring oechem.OEGraphMol is needed
    mol = oechem.OEGraphMol(mol)
    # add explicit hydrogens
    oechem.OEAddExplicitHydrogens(mol)
    mol_list.append(mol)
print(f"Loaded {len(mol_list)} molecules")

# write all molecules into separate mol2 files
folder_name = 'input_molecules'
if not os.path.isdir(folder_name):
    os.mkdir(folder_name)
os.chdir(folder_name)
if not os.path.isdir('xyz'):
    os.mkdir('xyz')
for idx, mol in enumerate(mol_list, 1):
    formula = oechem.OEMolecularFormula(mol)
    filename = f'{idx:03d}_{formula}.mol2'
    print(f'Writing {filename}')
    ofs.open(filename)
    oechem.OEWriteMolecule(ofs, mol)
    # write xyz files
    xyzfile = os.path.join('xyz', f'{idx:03d}_{formula}.xyz')
    ofs.open(xyzfile)
    oechem.OEWriteMolecule(ofs, mol)

