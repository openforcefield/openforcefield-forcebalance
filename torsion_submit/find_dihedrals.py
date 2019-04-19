#!/usr/bin/env python

import json
import collections
from forcebalance.molecule import Molecule, Elements

class DihedralSelector:
    def __init__(self, filename):
        self.m = Molecule(filename)
        self.build_neighbor_list()

    def build_neighbor_list(self):
        noa = self.m.na
        bonds = self.m.bonds
        neighbor_list = [[] for _ in range(noa)]
        for b1, b2 in bonds:
            neighbor_list[b1].append(b2)
            neighbor_list[b2].append(b1)
        for neighbors in neighbor_list:
            neighbors.sort()
        self.neighbor_list = neighbor_list

    def find_dihedrals(self, dihedral_filter="equiv"):
        """ find all dihedrals, filter out the ones with equivalent terminal atoms """
        noa = self.m.na
        # iterate over each bond, find distinct side atoms
        # indices i-j-k-l
        dihedral_list = []
        for j in range(noa):
            j_neighbors = self.neighbor_list[j]
            if len(j_neighbors) < 2: continue
            for k in j_neighbors:
                if k <= j: continue
                k_neighbors = self.neighbor_list[k]
                for i in j_neighbors:
                    if i == k: continue
                    for l in k_neighbors:
                        if l == j or l == i: continue
                        dihedral_list.append([i,j,k,l])
        print(f'Found total {len(dihedral_list)} distinct dihedrals')
        if dihedral_filter == "equiv":
            dihedral_list = self.filter_equivalent_terminals(dihedral_list)
        elif dihedral_filter == "heavy_no_ring":
            dihedral_list = self.filter_keep_4_heavy(dihedral_list)
            dihedral_list = self.filter_remove_ring(dihedral_list)
        return dihedral_list

    def filter_equivalent_terminals(self, dihedral_list):
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
        neighbor_list = self.neighbor_list
        noa = self.m.na
        equal_atom_idxs = {i:{i} for i in range(noa)}
        for i in range(noa):
            for j in range(i+1, noa):
                if elem_list[i] == elem_list[j]:
                    if len(neighbor_list[i]) == len(neighbor_list[j]) == 1:
                        if neighbor_list[i] == neighbor_list[j]:
                            equal_atom_idxs[i].add(j)
                            equal_atom_idxs[j].add(i)
        return equal_atom_idxs

    def filter_keep_4_heavy(self, dihedral_list):
        """ Filter dihedrals, keep only with 4 heavy atoms """
        print("Filter: only keep dihedrals formed by 4 heavy atoms")
        elem_list = self.m.elem
        filtered_dihedral_list = [d for d in dihedral_list if not any(elem_list[idx] == 'H' for idx in d)]
        print(f"Number Left: {len(dihedral_list)} => {len(filtered_dihedral_list)}")
        return filtered_dihedral_list

    def filter_remove_ring(self, dihedral_list):
        """ Filter dihedrals, remove any dihedral that's inside a ring """
        rings = self.find_atom_idxs_ring()
        print(f"Filter: Removing dihedrals that are containing in any of the rings {rings}")
        # build a dictionary stores the index of ring each atom belongs to
        d_rings = collections.defaultdict(set)
        for i_ring, ring in enumerate(rings):
            for atom_idx in ring:
                d_rings[atom_idx].add(i_ring)
        # go over dihedrals and check if any four belong to the same ring
        filtered_dihedral_list = []
        for i, j, k, l in dihedral_list:
            if d_rings[i] & d_rings[j] & d_rings[k] & d_rings[l]:
                print(f"Dihedral {i}-{j}-{k}-{l} skipped because in the same ring")
            else:
                filtered_dihedral_list.append([i,j,k,l])
        print(f"Number Left: {len(dihedral_list)} => {len(filtered_dihedral_list)}")
        return filtered_dihedral_list

    def find_atom_idxs_ring(self):
        """ Find rings in topology, return atom indices of each ring """
        distinct_rings = []
        paths = [[b0, b1] for b0, b1 in self.m.bonds]
        while len(paths) > 0:
            path = paths.pop()
            origin = path[0]
            last_idx = path[-1]
            for neighbor in self.neighbor_list[last_idx]:
                if neighbor in path:
                    if neighbor == origin and len(path) > 1:
                        # found a ring
                        if path[1] < last_idx:
                            # canonical ring index
                            # 1-2-3 will be kept, 1-3-2 will be skipped
                            distinct_rings.append(path)
                elif neighbor > origin:
                    # origin should be the smallest index in canonical ring
                    new_path = path + [neighbor]
                    paths.append(new_path)
        return distinct_rings

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
