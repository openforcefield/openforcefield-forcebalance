#!/usr/bin/env python

import json
import collections
from forcebalance.molecule import Molecule, Elements

class DihedralSelector:
    def __init__(self, filename):
        self.m = Molecule(filename)
        self.build_bond_list()

    def build_bond_list(self):
        noa = self.m.na
        bonds = self.m.bonds
        bond_list = [[] for _ in range(noa)]
        for b1, b2 in bonds:
            bond_list[b1].append(b2)
            bond_list[b2].append(b1)
        for neighbors in bond_list:
            neighbors.sort()
        self.bond_list = bond_list

    def find_dihedrals(self, filter=True):
        noa = self.m.na
        # iterate over each bond, find distinct side atoms
        # indices i-j-k-l
        dihedral_list = []
        for j in range(noa):
            j_neighbors = self.bond_list[j]
            if len(j_neighbors) < 2: continue
            for k in j_neighbors:
                if k <= j: continue
                k_neighbors = self.bond_list[k]
                for i in j_neighbors:
                    if i == k: continue
                    for l in k_neighbors:
                        if l == j or l == i: continue
                        dihedral_list.append([i,j,k,l])
        if filter is True:
            eq_idxs = self.find_equivalent_terminal_atom_idxs()
            print(f"Filtering based on equivalent terminal atoms: {eq_idxs}")
            dihedrals_set = set()
            for i,j,k,l in dihedral_list:
                skipped = False
                for eq_i in eq_idxs[i]:
                    for eq_l in eq_idxs[l]:
                        if (eq_i, j, k, eq_l) in dihedrals_set:
                            print(f"Filter: dihedral {i}-{j}-{k}-{l} skipped because equivalent to {eq_i}-{j}-{k}-{eq_l}")
                            skipped = True
                            break
                    if skipped: break
                if not skipped:
                    dihedrals_set.add((i,j,k,l))
            # resume the ordering of dihedrals
            dihedral_list = sorted(dihedrals_set, key=lambda d: (d[1], d[2], d[0], d[3]))
        return dihedral_list

    def find_equivalent_terminal_atom_idxs(self):
        elem_list = self.m.elem
        bond_list = self.bond_list
        noa = self.m.na
        equal_atom_idxs = {i:{i} for i in range(noa)}
        for i in range(noa):
            for j in range(i+1, noa):
                if elem_list[i] == elem_list[j]:
                    if len(bond_list[i]) == len(bond_list[j]) == 1:
                        if bond_list[i] == bond_list[j]:
                            equal_atom_idxs[i].add(j)
                            equal_atom_idxs[j].add(i)
        return equal_atom_idxs

    def write_dihedrals(self, dihedral_list, filename):
        with open(filename, 'w') as outfile:
            json.dump(dihedral_list, outfile, indent=2)

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help='Input mol2 file for a single molecule')
    parser.add_argument("-o", "--outfile", default='dihedral_list.json', help='Output json file containing definition of dihedrals')
    args = parser.parse_args()

    selector = DihedralSelector(args.infile)
    dihedral_list = selector.find_dihedrals()
    print(f"Found {len(dihedral_list)} dihedrals")
    for d in dihedral_list:
        print(d)
    selector.write_dihedrals(dihedral_list, args.outfile)

if __name__ == '__main__':
    main()
