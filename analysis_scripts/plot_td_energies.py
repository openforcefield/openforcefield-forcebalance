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

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def collect_td_targets_data(tmp_folder):
    """ Collect td targets QM vs MM data from tmp folder.
    Returns
    -------
    data: dict
        {
            'td_SMIRNOFF_Coverage_Torsion_Set_1_000_C3H8O3': {
                0: e_compare_data0, # iterations
                1: e_compare_data1,
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
        data[tgt_name] = {}
        tgt_folder = os.path.join(tmp_folder, tgt_name)
        iter_folders = [f for f in os.listdir(tgt_folder) if f.startswith('iter_') and os.path.isdir(os.path.join(tgt_folder, f))]
        for iter_name in iter_folders:
            iter_idx = int(iter_name[5:])
            data_file = os.path.join(tgt_folder, iter_name, 'EnergyCompare.txt')
            data[tgt_name][iter_idx] = np.loadtxt(data_file)
    return data

def plot_td_targets_data(data, folder_name='td_targets_plots', compare_first=True):
    ''' plot the target data, compare first and last iteration '''
    if os.path.exists(folder_name):
        shutil.rmtree(folder_name)
    os.mkdir(folder_name)
    print(f"Generating plots for {len(data)} targets in {folder_name}")
    for tgt_name, tgt_data in data.items():
        print(tgt_name)
        iter_list = sorted(tgt_data.keys())
        # create an new figure
        plt.Figure()
        last_iter_data = tgt_data[iter_list[-1]]
        # check qm - mm for last iter
        max_diff = np.max(np.abs(last_iter_data[:, 2]))
        if max_diff > 50:
            print(f"Warning, {tgt_name} iter {iter_list[-1]} max |qm-mm| > 50 kJ/mol")
        # plot qm
        plt.plot(last_iter_data[:,0], label='QM Relative Energies')
        # plot mm
        plt.plot(last_iter_data[:,1], label=f'MM Iter {iter_list[-1]}')
        # plot first iter
        if compare_first and len(iter_list) > 0:
            first_iter_data = tgt_data[iter_list[0]]
            plt.plot(first_iter_data[:, 1], label=f'MM Iter {iter_list[0]}')
        plt.legend()
        plt.xlabel('Frames')
        plt.ylabel('Relative Energies (kJ/mol)')
        plt.savefig(os.path.join(folder_name, tgt_name + '.pdf'))
        plt.close()




def main():
    import argparse
    parser = argparse.ArgumentParser('collect td targets data from fitting tmp folder and plot them')
    parser.add_argument('-f', '--tmp_folder', default='optimize.tmp')
    #parser.add_argument('-l', '--load_pickle', help='Load data directly from pickle file')
    #parser.add_argument('-i', '--iter', type=int, default=None, help='iteration number to read rmsd')
    args = parser.parse_args()

    data = collect_td_targets_data(args.tmp_folder)
    plot_td_targets_data(data)

if __name__ == '__main__':
    main()