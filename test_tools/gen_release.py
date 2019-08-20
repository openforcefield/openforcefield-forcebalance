#!/usr/bin/env python

import os
import shutil
import subprocess

this_file_folder = os.path.dirname(os.path.realpath(__file__))
analysis_script_folder = os.path.join(os.path.dirname(this_file_folder), 'analysis_scripts')
validate_script_folder = os.path.join(os.path.dirname(this_file_folder), 'validation_scripts')

# make a "release" folder
os.mkdir('release')
os.chdir('release')

# make a fb-fit folder
os.mkdir('fb-fit')
os.chdir('fb-fit')
print("release/fb-fit folder created")
print("Copying fitting related files")
# copy input/output/results into this folder
for fnm in ['optimize.in', 'optimize.out', 'smirnoff_parameter_assignments.json']:
    shutil.copyfile(os.path.join('../..', fnm), fnm)
for dnm in ['forcefield', 'targets', 'result']:
    shutil.copytree(os.path.join('../..', dnm), dnm)
os.chdir('..')

# copy the result force field file here
shutil.copyfile(os.path.join('fb-fit', 'result', 'optimize', 'param_valence.offxml'), 'result.offxml')


# make a "analysis" sub-folder
os.mkdir('analysis')
os.chdir('analysis')
print("release/analysis folder created")

# analysis script and running commands
analysis_scripts_cmd = {
    'visualize_fb_parameters.py': '../../optimize.out -x ../result.offxml',
    'plot_optgeo_each_smirks.py': '-x ../../forcefield/param_valence.offxml --new_xml ../result.offxml -f ../../optimize.tmp -t ../../targets -j ../../smirnoff_parameter_assignments.json',
    'plot_td_energies.py': '-f ../../optimize.tmp/ -t ../../targets/',
    'plot_vibfreq_rmsd.py': '-f ../../optimize.tmp',
}

def dualprint(msg, fp, **kwargs):
    """ print msg to both stdout and the file object """
    print(msg, **kwargs)
    print(msg, file=fp, flush=True, **kwargs)

with open('run_analysis.log', 'w') as fp:
    for script, arguments in analysis_scripts_cmd.items():
        command = ['python', os.path.join(analysis_script_folder, script)] + arguments.split()
        dualprint(f"\n\n-= Running script {script} =-", fp)
        dualprint(f"\n command:\n" + ' '.join(command), fp)
        print(f"\n stdout:\n", file=fp, flush=True)
        subprocess.run(command, check=True, stdout=fp, stderr=fp)

os.chdir('..')

print("All analysis results are generated successfully")

# make a "validate_tools" folder
os.mkdir('validate_tools')
os.chdir('validate_tools')
print("release/validate_tools folder created")

validate_scripts_cmd = {
    'prepare_validate_single_point.py': '-x ../result.offxml -t ../../targets',
}

with open('prep_validation_tools.log', 'w') as fp:
    for script, arguments in validate_scripts_cmd.items():
        command = ['python', os.path.join(validate_script_folder, script)] + arguments.split()
        dualprint(f"\n\n-= Running script {script} =-", fp)
        dualprint(f"\n command:\n" + ' '.join(command), fp)
        print(f"\n stdout:\n", file=fp, flush=True)
        subprocess.run(command, check=True, stdout=fp, stderr=fp)

with open('README.txt', 'w') as fnote:
    fnote.write("Validate tools are provided in this folder\n")
    fnote.write("To validate a test force field against the resulting force field in this release, go to the validate_single_point folder, then run\n\n")
    fnote.write('python run_validate_single_point.py test_forcefield.offxml\n')

print("All validation tools created successfully")

os.chdir('..')