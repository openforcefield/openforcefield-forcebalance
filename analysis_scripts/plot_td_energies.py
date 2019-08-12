#!/usr/bin/env python
"""
This script is used to analyze ForceBalance fitting resutls for optgeo targets.
Author: Yudong Qiu
"""

import os
import sys
import shutil
import numpy as np
import pickle
import json
from collections import Counter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def collect_td_targets_data(tmp_folder, targets_folder):
    """ Collect td targets QM vs MM data from tmp folder.
    Returns
    -------
    data: dict
        {
            'td_SMIRNOFF_Coverage_Torsion_Set_1_000_C3H8O3': {
                'metadata': {
                    'dihedrals': [[6, 10, 12, 11]],
                    'grid_spacing': [15],
                    'dihedral_ranges': None,
                    'energy_decrease_thresh': None,
                    'energy_upper_limit': 0.05,
                    'dataset_name': 'OpenFF Group1 Torsions',
                    'entry_label': 'c1c[cH:1][c:2](cc1)[CH2:3][c:4]2ccccc2',
                    'canonical_smiles': 'c1ccc(cc1)Cc2ccccc2',
                    'torsion_grid_ids': [[-165], [-150], [-135], [-120], ...],
                    'smirks': ['[*:1]~[#6X3:2]-[#6X4:3]-[*:4]'],
                    'smirks_ids': ['t17'],
                },
                'iterdata': {
                    0: e_compare_data0,
                    1: e_compare_data1,
                }

            }
        }
        where e_compare_data is a 2D numpy array [[eqm, emm, diff, weight], ..]
    """
    td_target_folders = [f for f in os.listdir(tmp_folder) if f.startswith('td_') and os.path.isdir(os.path.join(tmp_folder, f))]
    td_target_folders.sort()
    print(f'Collecting td data from {len(td_target_folders)} target folders')
    data = {}
    for i, tgt_name in enumerate(td_target_folders):
        print(i, tgt_name)
        # load metadata from targets folder
        with open(os.path.join(targets_folder, tgt_name, 'metadata.json')) as jsonfile:
            metadata = json.load(jsonfile)
        data[tgt_name] = {'metadata': metadata, 'iterdata': {}}
        tgt_folder = os.path.join(tmp_folder, tgt_name)
        iter_folders = [f for f in os.listdir(tgt_folder) if f.startswith('iter_') and os.path.isdir(os.path.join(tgt_folder, f))]
        for iter_name in iter_folders:
            iter_idx = int(iter_name[5:])
            data_file = os.path.join(tgt_folder, iter_name, 'EnergyCompare.txt')
            data[tgt_name]['iterdata'][iter_idx] = np.loadtxt(data_file, ndmin=2)
    return data

def plot_td_targets_data(data, folder_name='td_targets_plots', compare_first=True, iteration=None):
    ''' plot the target data, compare first and last iteration '''
    # aggregate metadata for counts
    smirks_counter = Counter(sid for tgt_data in data.values() for sid in tgt_data['metadata']['smirks_ids'])
    # create new folder
    if os.path.exists(folder_name):
        shutil.rmtree(folder_name)
    os.mkdir(folder_name)
    print(f"Generating plots for {len(data)} targets in {folder_name}")
    for tgt_name, tgt_data in data.items():
        print(tgt_name)
        metadata = tgt_data['metadata']
        iterdata = tgt_data['iterdata']
        iter_list = sorted(iterdata.keys())
        # create an new figure
        plt.Figure()
        last_iter = iteration if iteration != None else iter_list[-1]
        last_iter_data = iterdata[last_iter]
        # check qm - mm for last iter
        max_diff = np.max(np.abs(last_iter_data[:, 2]))
        if max_diff > 50:
            print(f"Warning, {tgt_name} iter {last_iter} max |qm-mm| > 50 kJ/mol")
        # plot qm
        plt.plot(last_iter_data[:,0], label='QM Relative Energies', color='C0')
        # plot mm
        plt.plot(last_iter_data[:,1], label=f'MM Iter {last_iter}', color='C2')
        # plot first iter
        if compare_first and len(iter_list) > 0:
            first_iter_data = iterdata[iter_list[0]]
            plt.plot(first_iter_data[:, 1], label=f'MM Iter {iter_list[0]}', color='C1')
        plt.legend()
        # use strings as x axis to support >1D torsion scans
        torsion_labels = [str(gid)[1:-1] for gid in metadata['torsion_grid_ids']]
        tick_x_locs = list(range(len(last_iter_data)))
        tick_stripe = max(int(len(last_iter_data)/10), 1) # reduce tick density
        plt.xticks(ticks=tick_x_locs[::tick_stripe], labels=torsion_labels[::tick_stripe])
        # print metadata as footnote
        footnotes = {
            "Dataset Name": metadata['dataset_name'],
            "Entry Label": metadata['entry_label'],
            "Canonical SMILES": metadata['canonical_smiles'],
            "Torsion Atom Indices": ', '.join(map(str, metadata['dihedrals'])),
            "Torsion SMIRKs": ', '.join(metadata['smirks']),
            "Torsion SMIRKs ID": ', '.join(metadata['smirks_ids']),
            "SMIRKs Total Count": ', '.join(str(smirks_counter[sid]) for sid in metadata['smirks_ids']),
        }
        footnote_lines = [f'{key:<20s} {value:>59s}' for key, value in footnotes.items()]
        plt.xlabel('Torsion Angles\n\n' + '\n'.join(footnote_lines), fontdict={'family':'monospace', 'size': 8})
        # plt.xlabel('Torsion Angles')
        plt.ylabel('Relative Energies (kcal/mol)')
        plt.tight_layout()
        filename = os.path.join(folder_name, tgt_name + '.pdf')
        plt.savefig(filename)
        plt.close()
        # copy the figure to subfolders of the target torsions SMIRKs
        for sid in metadata['smirks_ids']:
            subfolder = os.path.join(folder_name, sid)
            if not os.path.exists(subfolder):
                os.mkdir(subfolder)
            # move file if it's the last one
            if sid == metadata['smirks_ids'][-1]:
                shutil.move(filename, os.path.join(subfolder, tgt_name + '.pdf'))
            else:
                shutil.copyfile(filename, os.path.join(subfolder, tgt_name + '.pdf'))





def main():
    import argparse
    parser = argparse.ArgumentParser('collect td targets data from fitting tmp folder and plot them')
    parser.add_argument('-f', '--tmp_folder', default='optimize.tmp')
    parser.add_argument('-t', '--targets_folder', default='targets')
    parser.add_argument('-l', '--load_pickle', help='Load data directly from pickle file')
    parser.add_argument('-i', '--iteration', type=int, default=None, help='iteration number to read rmsd')
    args = parser.parse_args()

    if args.load_pickle:
        with open(args.load_pickle, 'rb') as pfile:
            data = pickle.load(pfile)
    else:
        data = collect_td_targets_data(args.tmp_folder, args.targets_folder)
        # save data as pickle file
        with open('td_plot_data.pickle', 'wb') as pfile:
            pickle.dump(data, pfile)
    plot_td_targets_data(data, iteration=args.iteration)

if __name__ == '__main__':
    main()