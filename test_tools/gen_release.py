#!/usr/bin/env python

import os
import shutil
import subprocess

this_file_folder = os.path.dirname(os.path.realpath(__file__))
analysis_script_folder = os.path.join(os.path.dirname(this_file_folder), 'analysis_scripts')

# make a "release" folder
os.mkdir('release')
os.chdir('release')

# make a fb-fit folder
os.mkdir('fb-fit')
os.chdir('fb-fit')
# copy input/output/results into this folder
for fnm in ['optimize.in', 'optimize.out', 'smirnoff_parameter_assignments.json']:
    shutil.copyfile(os.path.join('../..', fnm), fnm)
for dnm in ['forcefield', 'targets', 'result']:
    shutil.copytree(os.path.join('../..', dnm), dnm)
os.chdir('..')
print("release/fb-fit folder created successfully")

# make a "analysis" sub-folder
os.mkdir('analysis')
os.chdir('analysis')

# analysis script and running commands
analysis_scripts_cmd = {
    'visualize_fb_parameters.py': '../../optimize.out -x ../../forcefield/param_valence.offxml',
    'plot_optgeo_each_smirks.py': '-x ../../forcefield/param_valence.offxml --new_xml ../../result/optimize/param_valence.offxml -f ../../optimize.tmp -t ../../targets -j ../../smirnoff_parameter_assignments.json',
    'plot_td_energies.py': '-f ../../optimize.tmp/ -t ../../targets/',
    'plot_vibfreq_rmsd.py': '-f ../../optimize.tmp',
}

def dualprint(msg, fp, **kwargs):
    """ print msg to both stdout and the file object """
    print(msg, **kwargs)
    print(msg, file=fp, **kwargs)

with open('run_analysis.log', 'w') as fp:
    for script, arguments in analysis_scripts_cmd.items():
        command = ['python', os.path.join(analysis_script_folder, script)] + arguments.split()
        dualprint(f"\n\n-= Running script {script} =-", fp)
        dualprint(f"\n command:\n" + ' '.join(command), fp)
        print(f"\n stdout:\n", file=fp)
        subprocess.run(command, check=True, stdout=fp)

os.chdir('..')

print("release/analysis folder created successfully")