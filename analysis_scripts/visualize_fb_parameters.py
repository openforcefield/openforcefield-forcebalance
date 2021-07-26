#!/usr/bin/env python

import numpy as np
import collections
import matplotlib
matplotlib.use('Agg')
matplotlib.rc('ytick', labelsize=6)
import matplotlib.pyplot as plt
import sys

def read_fb_params(filename):
    readlines = open(filename).readlines()
    param_list = collections.OrderedDict()
    # count the number of steps found
    n_steps = 0
    # find the initial paramters
    found_initial = False
    for line in readlines:
        if 'Starting parameter indices, physical values and IDs' in line:
            found_initial = True
            n_steps += 1
        if found_initial:
            ls = line.split()
            if len(ls) == 6:
                idx, value = ls[0], float(ls[2])
                ptype = '/'.join(ls[5].split('/')[:-1])
                #ForceBalance output file lists pytpes differently between GROMACS and OpenMM
                #Uncomment the below line to work with GROMACS output
                #ptype = ls[5].split(':')[0]
                name = idx + ': ' + ls[5].split('/')[-1]
                if ptype not in param_list:
                    param_list[ptype] = dict()
                param_list[ptype][name] = [value]
            elif line.startswith("------"):
                break
    # read the prior width
    found_prior_width = False
    for line in readlines:
        if 'Rescaling Factors by Type' in line:
            found_prior_width = True
        if found_prior_width:
            ls = line.split()
            if len(ls) == 3:
                prior_type, value = ls[0], float(ls[2])
                for ptype in param_list:
                    if ptype.startswith(prior_type):
                        #print("Prior: setting %s to %f"%(ptype, value))
                        param_list[ptype]['PriorWidth'] = [value]
            elif line.startswith("------"):
                break
    #print param_list['Proper/k1']
    # read the paramters from each step
    found_step = False
    for line in readlines:
        if 'Physical Parameters (Current + Step = Next)' in line:
            found_step = True
            n_steps += 1
        if found_step:
            ls = line.split()
            if len(ls) == 10:
                idx, value = ls[0], float(ls[6])
                #ForceBalance output file lists pytpes differently between GROMACS and OpenMM
                #Uncomment the below line to work with GROMACS output
                ptype = '/'.join(ls[9].split('/')[:-1])
                #ptype = ls[9].split(':')[0]
                name = idx + ': ' + ls[9].split('/')[-1]
                if ptype not in param_list:
                    param_list[ptype] = dict()
                if name not in param_list[ptype]:
                    param_list[ptype][name] = []
                param_list[ptype][name].append(value)
            if line.startswith("------"):
                found_step = False
    # check we have the same number of steps for each parameter
    for ptype in param_list:
        for name in param_list[ptype]:
            if name != 'PriorWidth' and len(param_list[ptype][name]) != n_steps:
                raise RuntimeError("Inconsistent number of steps for %s/%s"%(ptype, name))
    return param_list


def plot_paramters(param_list):
    for ptype in param_list:
        # sort the parameters by their changes
        names, values = [], []
        prior_value = None
        for name, value in sorted(param_list[ptype].items(), key=lambda x:abs(x[1][-1] - x[1][0])):
            if name != 'PriorWidth':
                names.append(name)
                values.append(value)
            else:
                prior_value = value
        initial_values = [v[0] for v in values]
        value_changes = [v[-1] - v[0] for v in values]
        # add the prior width bar
        if prior_value is not None:
            names = ["Prior Width"] + names
            initial_values = [0] + initial_values
            value_changes = [0] + value_changes
        initial_values = np.array(initial_values)
        value_changes = np.array(value_changes)
        # the position of each bar
        y_pos = np.arange(len(names))
        # adjust the size of the figure
        plt.figure(figsize=(8.5,len(names)*0.12+0.8))
        # linewidth
        lw = None
        # plot the initial parameters
        final_values = initial_values + value_changes
        xmin = min(initial_values.min(), final_values.min())
        xmax = max(initial_values.max(), final_values.max())
        head_length = 0.01*(xmax-xmin)
        padding = (xmax - xmin) * 0.01
        plt.xlim(xmin, xmax+padding)
        if prior_value is not None:
            initial_values[0] = prior_value[0]
        plt.scatter(initial_values, y_pos, marker='o', s=80, facecolors='none', edgecolors='grey')
        plt.yticks(y_pos, names)
        # plot the changes in the final parameters
        increase_idxs = np.nonzero(value_changes >=0)[0]
        decrease_idxs = np.nonzero(value_changes <0)[0]
        for i in increase_idxs:
            if abs(value_changes[i]) > head_length:
                plt.arrow(initial_values[i], y_pos[i], value_changes[i], 0.0, head_width=0.4, head_length=head_length, length_includes_head=True, width=0.05, color='C3')
        for i in decrease_idxs:
            if abs(value_changes[i]) > head_length:
                plt.arrow(initial_values[i], y_pos[i], value_changes[i], 0.0, head_width=0.4, head_length=head_length, length_includes_head=True, width=0.05, color='C3')
        ## plot the prior width
        plt.title(ptype)
        # adjust the y range
        plt.ylim(y_pos[-1]+1, y_pos[0]-1)
        # adjust the x range
        plt.tight_layout()
        filename = ptype.replace('/', '_') + '.pdf'
        plt.savefig(filename)
        plt.close()

def main():
    import argparse, os, shutil, subprocess
    parser = argparse.ArgumentParser(description='Read a ForceBalance output file, bar plot the difference between initial and final parameters.')
    parser.add_argument('fbout', help='ForceBalance output file')
    parser.add_argument('-o', '--outfolder', default='param_change', help='Folder name to save each plot as a separate file')
    args = parser.parse_args()


    param_list = read_fb_params(args.fbout)
    if os.path.exists(args.outfolder):
        shutil.rmtree(args.outfolder)
    os.mkdir(args.outfolder)
    os.chdir(args.outfolder)
    plot_paramters(param_list)
    subprocess.call('pdfunite *.pdf all.pdf', shell=True)

if __name__ == '__main__':
    main()
