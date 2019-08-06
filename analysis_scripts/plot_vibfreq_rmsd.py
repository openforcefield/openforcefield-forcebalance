#!/usr/bin/env python
"""
This script is used to analyze ForceBalance fitting resutls for optgeo targets.
Author: Yudong Qiu
"""

import os
import shutil
import numpy as np
import pickle

import matplotlib
matplotlib.use('Agg')
matplotlib.rc('ytick', labelsize=6)
import matplotlib.pyplot as plt

def collect_vibfreq_data(tmp_folder):
    """ Collect td targets QM vs MM data from tmp folder.
    Returns
    -------
    data: dict
        {
            'td_SMIRNOFF_Coverage_Torsion_Set_1_000_C3H8O3': {
                'iterdata': {
                    0: freq_compare_data0,
                    1: freq_compare_data0,
                }

            }
        }
        where e_compare_data is a 2D numpy array [[freq_qm, freq_mm], ..]
    """
    td_target_folders = [f for f in os.listdir(tmp_folder) if f.startswith('vibfreq_') and os.path.isdir(os.path.join(tmp_folder, f))]
    td_target_folders.sort()
    print(f'Collecting td data from {len(td_target_folders)} target folders')
    data = {}
    for i, tgt_name in enumerate(td_target_folders):
        print(i, tgt_name)
        data[tgt_name] = {'iterdata': {}}
        tgt_folder = os.path.join(tmp_folder, tgt_name)
        iter_folders = [f for f in os.listdir(tgt_folder) if f.startswith('iter_') and os.path.isdir(os.path.join(tgt_folder, f))]
        for iter_name in iter_folders:
            iter_idx = int(iter_name[5:])
            data_file = os.path.join(tgt_folder, iter_name, 'indicate.log')
            data[tgt_name]['iterdata'][iter_idx] = read_indicate(data_file)
    return data

def read_indicate(fnm):
    """ Read reference and calculated frequencies """
    res = []
    with open(fnm) as f:
        for line in f:
            if line[0] != '#':
                ls = line.split()
                if len(ls) == 5:
                    ref = float(ls[1])
                    calc = float(ls[2])
                    res.append([ref, calc])
    return np.array(res, dtype=float)

def get_rmsd_maxdiff(freq_compare_data):
    """ Get the RMSD and maxdiff value from the freq_compare_data """
    diff = freq_compare_data[:, 0] - freq_compare_data[:, 1]
    sq_diff = np.square(diff)
    max_diff = np.sqrt(np.max(sq_diff))
    rmsd = np.sqrt(np.sum(sq_diff) / len(sq_diff))
    return rmsd, max_diff

def plot_vibfreq_data(data, folder_name='vibfreq_targets_plots', iteration=None):
    ''' plot the target data, compare first and last iteration '''
    # determine the last iteration to plot
    target_names = list(data)
    iter_list = sorted(data[target_names[0]]['iterdata'])
    last_iter = iteration if iteration != None else iter_list[-1]
    # compute vibfreq rmsd and max_diff
    tgt_rmsd_list = []
    tgt_max_diff_list = []
    tgt_rmsd_list_iter0 = []
    tgt_max_diff_list_iter0 = []
    for tgt_name, tgt_data in data.items():
        rmsd, maxdiff = get_rmsd_maxdiff(tgt_data['iterdata'][last_iter])
        tgt_rmsd_list.append(rmsd)
        tgt_max_diff_list.append(maxdiff)
        # compute for the first iter
        rmsd_iter0, maxdiff_iter0 = get_rmsd_maxdiff(tgt_data['iterdata'][0])
        tgt_rmsd_list_iter0.append(rmsd_iter0)
        tgt_max_diff_list_iter0.append(maxdiff_iter0)

    # create new folder
    if os.path.exists(folder_name):
        shutil.rmtree(folder_name)
    os.mkdir(folder_name)
    print(f"Generating rmsd and maxdiff plots for {len(data)} targets in folder {folder_name}")
    # make the rmsd plot
    title = "QM vs MM RMSD of vibration frequencies"
    xlabel = r'RMSD of Vibrational Frequencies ( cm $^{-1}$ )'
    fnm = os.path.join(folder_name, 'vib_freq_rmsd.pdf')
    gen_plot(fnm, target_names, tgt_rmsd_list_iter0, tgt_rmsd_list, last_iter, title, xlabel)
    # make the maxdiff plot
    title = "Max QM vs MM difference of vibration frequencies"
    xlabel = r'Max QM vs MM Differences of Vibrational Frequencies ( cm $^{-1}$ )'
    fnm = os.path.join(folder_name, 'vib_freq_maxdiff.pdf')
    gen_plot(fnm, target_names, tgt_max_diff_list_iter0, tgt_max_diff_list, last_iter, title, xlabel)

def gen_plot(filename, target_names, data_iter0, data_last_iter, last_iter, title, xlabel):
    # create a new figure
    fig = plt.figure(figsize=(8.5,len(target_names)*0.12+0.8))
    y_pos = np.arange(len(target_names))[::-1]
    # linewidth
    lw = None
    # plot the initial rmsd
    plt.barh(y_pos, data_iter0, height=0.4, color='C1', tick_label=target_names, align='edge', linewidth=lw, label='iter 0')
    # plot the final rmsd
    plt.barh(y_pos-0.2, data_last_iter, height=0.4, color='C2', align='center', linewidth=lw, label=f'iter {last_iter}')
    # adjust the y range
    plt.ylim(y_pos.min()-1, y_pos.max()+1)
    plt.legend()
    plt.title(title)
    plt.xlabel(xlabel)
    # adjust the y tick labels to align left
    plt.yticks(ha='left')
    ax = plt.gca()
    ax.yaxis.set_tick_params(pad=180)
    # save plot
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

def main():
    import argparse
    parser = argparse.ArgumentParser('collect vibfreq targets data from fitting tmp folder and plot them')
    parser.add_argument('-f', '--tmp_folder', default='optimize.tmp')
    parser.add_argument('-l', '--load_pickle', help='Load data directly from pickle file')
    parser.add_argument('-i', '--iteration', type=int, default=None, help='iteration number to read rmsd')
    args = parser.parse_args()

    if args.load_pickle:
        with open(args.load_pickle, 'rb') as pfile:
            data = pickle.load(pfile)
    else:
        data = collect_vibfreq_data(args.tmp_folder)
        # save data as pickle file
        with open('vibfreq_plot_data.pickle', 'wb') as pfile:
            pickle.dump(data, pfile)
    plot_vibfreq_data(data, iteration=args.iteration)

if __name__ == '__main__':
    main()