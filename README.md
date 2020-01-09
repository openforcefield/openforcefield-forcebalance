# OpenFF ForceBalance Fitting Tutorial

**Note 1**: Please read the [setup guide](setup_guide.md) for the installation of ForceBalance and other packages. 
**Note 2**: The fitting results are posted as [release packages](https://github.com/openforcefield/openforcefield-forcebalance/releases) in this repository. You can find `release_XX.tar.gz` which contains fitting results and all codes used in the fitting procedure for the release.

## Overview
Fitting valence parameters to QM data can be performed in four steps. 

(1) QM data generation and submission. Quantum chemical calculations, which are the most expensive part of the fitting procedure, are performed using the [MolSSI QCFractal](https://qcarchive.molssi.org/) distributed quantum chemistry engine.

(2) Pulling QM data from QCArchive and processing into ForceBalance native formats. The data is stored as `DataSet` objects on the public MolSSI QCArchive server and can be pulled using the `qcportal` Python API. In this step, we pull the data, filter and format them into ForceBalance fitting targets. Typical filters include “toolkit capability”, “topology consistency”, and “avoid hydrogen bonds”.

(3) Optimization of valence parameters by fitting to QM data. The main fitting procedure is carried out by the open-source software package [ForceBalance](https://github.com/leeping/forcebalance), which is maintained by the [Wang research group](https://www.lpwchem.org). In this step we prepare the input files for the ForceBalance command line application, then run ForceBalance to carry out the optimization.

(4) Analysis of the optimized force field. After the fitting is finished, the output file and temporary folder contain many details; to get insights into these, we perform several analyses that aggregate the fitting data and produce various plots and tables. 

## 4 Steps of Valance Parameter Fitting
### 1. QM data submission 
Quantum chemical calculations are performed on a distributed set of high-performance computing clusters using the MolSSI QCFractal distributed quantum chemistry engine, and stored as `DataSet` objects on the public MolSSI QCArchive server. 

For the release-1 "Parsley" parameter optimization, two sets of molecules were used, and each leads to three types of QM data being generated: 1) Optimized geometries; 2) Vibrational frequencies; 3) Torsion profiles.  Thus there are 6 QM data sets in total. 

**Table 1**. Collection types and names of datasets of 6 QM data sets specified on QCArchive server.
<center>

|                         |    Collection type    |    Name for "Roche set"   |      Name for "Coverage set"      |
|-------------------------|:---------------------:|:-------------------------:|:---------------------------------:|
| Optimized geometries    |  OptimizationDataset  | OpenFF Optimization Set 1 |  SMIRNOFF Coverage Set 1          |
| Vibrational frequencies |  Dataset              | OpenFF Optimization Set 1 |  SMIRNOFF Coverage Set 1          |
| Torsion profiles        |  TorsionDriveDataset  | OpenFF Group1 Torsions    |  SMIRNOFF Coverage Torsion Set 1  |  
</center>

The [blog post](https://openforcefield.org/news/introducing-openforcefield-1.0/#fitting-parsley-to-quantum-chemical-data) provides a nice overview of the background behind the selection and generation of datasets.
For details about how to generate and submit data, please check https://github.com/openforcefield/qca-dataset-submission.


### 2. Pulling QM data from QCArchive and processing into ForceBalance native formats for parameter fitting.
**Note: The result of this step  was aggregated into the fb-fit/targets/ folder of the release package. Thus if you want to reproduce the fitting without changing the targets, please skip to the next step.**

The results of quantum chemical calculations are stored as `DataSet` objects on the public MolSSI QCArchive and can be pulled from the server.

To list the available QM data sets on QCArchive server, run the following in a Jupyter notebook or Python console:
```
import qcportal as ptl
client = ptl.FractalClient('https://api.qcarchive.molssi.org:443/')
client.list_collections()
```
The list of all available collections currently on the server will be printed to the console. 

#### 2.1 Preparation
Open a terminal and go to an empty folder. 
From now on, all related files will be placed into this folder.

1) Clone the following repository for a set of scripts to convert the QCArchive data:

    ```
    git clone https://github.com/openforcefield/openforcefield-forcebalance.git
    ````

2) Create a new folder to carry out the remaining steps:

    ```
    mkdir make-target
    cd make-target
    ```

3) To ensure the targets generated are compatible with the toolkit using this force field file, we copy over the force field `.offxml` file that is going to be used in parameter fitting. One such file can be found in the release package `fb-fit/forcefield/param_valence.offxml`. You can also download the input force field file we used in the release-1 fitting:

    ```
    wget https://github.com/openforcefield/openforcefield-forcebalance/releases/download/v1.0.0-RC2/input_ff.offxml
    ```

#### 2.2. Prepare optimized geometry fitting targets 
In release-1 we used two OptimizationDatasets, namely “OpenFF Optimization Set 1” and “SMIRNOFF Coverage Set 1”. To pull “OpenFF Optimization Set 1” and generate a list of ForceBalance targets:
```
mkdir optgeo-opt-set1; cd optgeo-opt-set1
python ~/codes/openforcefield-forcebalance/optimized_geo/make_fb_optgeo_target.py “OpenFF Optimization Set 1” -t ../input_ff.offxml | tee run.log
```
This will 1) pull the data set from QCArchive; 2) filter the data to exclude cases where the bonding pattern changes during QM optimization and cases where the data introduces errors for the openforcefield toolkit. The molecules which are filtered out can be found in error_mol2s in the target folder generated. and 3) format the data into ForceBalance targets.

The progress will be printed on the screen and saved in the log file `run.log`. The resulting files and folders are saved in the `targets/` folder. Each subfolder is a ForceBalance target that we will use later in the fitting. A file named “targets/target.XXXX.in” is also generated that will be used when preparing the ForceBalance input file.

Repeat the same commands for the “Coverage set”:
```
cd ..; mkdir optgeo-coverage-set1; cd optgeo-coverage-set1
python ../../openforcefield-forcebalance/optimized_geo/make_fb_optgeo_target.py “SMIRNOFF Coverage Set 1” -t ../input_ff.offxml | tee make_target.log
```
Notice that the only argument changed is the name of the dataset “SMIRNOFF Coverage Set 1”.

#### 2.3. Prepare vibrational frequency targets
To pull the hessian data from the DataSet “OpenFF Optimization Set 1” and generate vibrational frequency fitting targets:
```
cd ..; mkdir vibfrq-opt-set1; cd vibfrq-opt-set1
python ../../openforcefield-forcebalance/vib_freq_target/make_vib_frq_target.py “OpenFF Optimization Set 1” -t ../input_ff.offxml | tee run.log
```
This will 1) download Hessian data for each optimized geometry; 2) pick the lowest-energy conformer of each molecule; 3) apply the toolkit-compatibility and topology-check filters; 4) perform normal mode analysis to get the vibrational frequencies; and 5) write the data into `targets/` folder for ForceBalance.

To repeat the process for “SMIRNOFF Coverage Set 1”:
```
cd ..; mkdir vibfrq-coverage-set1; cd vibfrq-coverage-set1
python ../../openforcefield-forcebalance/vib_freq_target/make_vib_frq_target.py “SMIRNOFF Coverage Set 1” -t ../input_ff.offxml | tee runmake_target.log
```

#### 2.4. Prepare the torsion profile fitting targets 
To generate the fitting targets for a TorsionDriveDataset “OpenFF Optimization Set 1”:
```
cd ..; mkdir td-opt-set1; cd td-opt-set1
python ../../openforcefield-forcebalance/torsion_target/make_torsion_target_new.py “OpenFF Group1 Torsions” -t ../input_ff.offxml | tee run.log
```
This script pulls torsion scan trajectories from the server, filters out trajectories that contain any frame with hydrogen bonds, and formats into ForceBalance torsion profiles targets, while saving the useful metadata about the torsiondrive records as metadata.json in each target folder. 

Repeating the process for “SMIRNOFF Coverage Set 1”:
```
cd ..;mkdir td-coverage-set1;cd td-coverage-set1
python ../../openforcefield-forcebalance/torsion_target/make_torsion_target_new.py “SMIRNOFF Coverage Torsion Set 1” -t ../input_ff.offxml | tee run.log
```

### 3. Run Fitting with ForceBalance
ForceBalance is a free software designed for carrying out force field optimizations using a systematic and reproducible procedure. 

Note: The result of steps 3.1 and 3.2 are provided in the release package.

#### 3.1. Create fb-fit in the same location as the `make-target/` folder from last step. We are going to run ForceBalance in here:
```
mkdir fb-fit
cd fb-fit
```
#### 3.2. Prepare three input components for ForceBalance 
Notice that all needed files are already provided in the release package under fb-fit/ folder. The following steps can be used to re-create these files.

#### 3.2.1 Input force field 

The directory named `forcefield` contains the force field files that you are optimizing. Parameters to be optimized are specified by `parameterize` tags. 
For example, if we want to fit the force constant and the equilibrium bond length of the bond term with string ID `b1`, the attribute `parameterize=”k,length”` should be added as shown below.
```
<Bond smirks="[#6X4:1]-[#6X4:2]" length="1.526 * angstrom" k="620.0 * angstrom**-2 * mole**-1 * kilocalorie" id="b1" parameterize="k,length"/>
```

You can get the directory directly from the release package and use it for the calculation, or you can manually modify the input force field files as needed. 
When doing so, it is preferred to use the openforcefield toolkit for the modification of the file.

#### 3.2.2 Fitting targets

The directory named `targets` contains reference data sets as well as input files for simulating that data using the force field. 
Each subdirectory in `targets` corresponds to a single "ForceBalance target" or data set, and its contents depends on the specific target type. 
Since the targets from the previous step can be directly used in ForceBalance fitting, we can copy and merge all the prepared files into one folder.
```
cp -r ../make-target/*/targets .
```

#### 3.2.3 Input file

The ForceBalance input file contains a "global options" section, starting with `$options`, that specifies global properties of the optimization.
There are also "target options" section(s), starting with `$target`, which specify the settings of one or more ForceBalance targets.

An example input file is provided as `optimize.in` in the release package under the `fb-fit/` folder. 

Global options: 

`jobtype`, specifies the calculation type. To run force field optimization, set jobtype to `optimize`. To run single-point calculation, set jobtype to `single`. 

The provided `optimize.in` contains two global options `wq_port`, and `asynchronous` that specify the use of Work Queue for distributed target evaluations. 
Please refer to the setup guide if you are using Work Queue for the first time. 
If you intend to run small optimization jobs on your local machine, please remove these two options. 

The settings for fitting targets are given as `$target` blocks, one for each target. 
The provided `optimize.in` contains a long list of such blocks. 
To reproduce these input blocks, you can clean upall of them from `optimize.in`, then append the updated input options by:
```
cat targets/target*.in >> optimize.in
```

#### 3.3. Run ForceBalance 

The ForceBalance optimizer is executed on the command line as:
```
ForceBalance optimize.in
```
The executable is intended to be run in the foreground, allowing you to view outputs in real time. 
The output will be printed on the screen and also saved in the output file `optimize.out`. 
Although it is possible to run in the background, the recommended workflow is to run ForceBalance in an open terminal window or a virtual console, e.g. provided by the GNU Screen or tmux software packages.

All fitting related temporary files containing the details are saved in the directory `optimize.tmp/`. 
This folder is too large to be included in the release package. 
After the fitting is finished, the resulting forcefield file is saved in the `result/` folder.

### 4. Analysis of ForceBalance optimization result
Several scripts for visualization of fitting results can be found in `analysis_scripts/` directory. They help to get insights on ForceBalance output data. 

Note: The results of this step is provided in the release package in the `analysis/` folder.

#### 4.1 Create a new folder for saving the analysis results, in the same location as the fb-fit/ folder.
```
mkdir analysis
cd analysis
```

#### 4.2. Visualize the change of parameters in fitting
In this step, we will read the ForceBalance output file to generate bar plots of parameter changes, with green and red colors indicating an increase/decrease in the parameter value respectively. 
```
python ../openforcefield-forcebalance/analysis_scripts/visualize_fb_parameters.py ../fb-fit/optimize.out -x ../fb-fit/forcefield/param_valence.offxml 
```

The result of this script is saved in a folder `param_change/`, with each type of parameter in one pdf file. 
The `all.pdf` combines the parameter changes for all parameter types into a single file.

#### 4.3. Visualizing optimization results for optimized geometries
To visualize the effect of fitting on optimized geometries, we run a script to read the QM and MM internal coordinate values for every molecule from the `optimize.tmp` folder.
The script aggregates the internal coordinate values according to their matching SMIRKS patterns from the force field file and draws scatter plots with initial and final values of the equilibrium parameter as vertical lines. 
Each dot corresponds to one internal coordinate in an optimized structure.
```
python ../openforcefield-forcebalance/analysis_scripts/plot_optgeo_each_smirks.py -x ../fb-fit/forcefield/param_valence.offxml --new_xml../fb-fit/result/optimize/param_valence.offxml -f ../fb-fit/optimize.tmp -t ../fb-fit/targets-j ../fb-fit/smirnoff_parameter_assignments.json
```
This script will generate a folder `optgeo_scatter_plots/`, that contains subfolders like `Bonds/`. 
The PDF files in each folder correspond to the SMIRKS id, such as `b1.pdf`. 
This script will also save a pickle file `optgeo_analysis_data.p` on disk for quickly generating the plots without re-reading the data from `optimize.tmp`.

#### 4.4. Visualizing optimization results for vibrational frequencies
To visualize the effect of fitting on vibrational frequencies, a script was implemented to plot the bar chart of RMSD of QM and MM vibrational frequencies from the zeroth optimization cycle (before any parameter changes) and the final optimization cycle. 
```
python ../openforcefield-forcebalance/analysis_scripts/plot_vibfreq_rmsd.py -f ../fb-fit/optimize.tmp
```

This script will generate a folder `vibfreq_targets_plots/` that contains two files, the `vib_freq_rmsd.pdf` and `vib_freq_maxdiff.pdf`. 
A pickle file `vibfreq_plot_data.pickle` is also saved for re-using the loaded data when making new plots.

#### 4.5. Visualizing optimization results for torsion profiles
We use a script to load data from `optimize.tmp/` directory and generate QM and MM energy profiles as a function of the the torsion angle, including the MM energies from both the zeroth cycle (before any parameter changes) and the final optimization cycle for comparison. 
A table of metadata is also included in each plot. 
Plot files are organized into subfolders by the ID of the matching SMIRKS pattern.
```
python ../openforcefield-forcebalance/analysis_scripts/plot_td_energies.py -f ../fb-fit/optimize.tmp -t ../fb-fit/targets
```

This script will generate a folder `td_targets_plots/` that contains subfolders corresponding to each SMIRKs torsion term. 
Each PDF file has the same name as the corresponding torsion profile target such as `td_OpenFF_Group1_Torsions_193_C11H18N2.pdf`.
The log of running all above analysis scripts is saved in the release package as `analysis/run_analysis.log`.

## End
Congratulations! Now you have completed the full fitting tutorial. 
Please make sure you read the [blog post](https://openforcefield.org/news/introducing-openforcefield-1.0/#fitting-parsley-to-quantum-chemical-data) for the theory and ideas behind this optimization.

For questions about this tutorial, details on the fitting releases, or ForceBalance, please contact Hyesu Jang, Yudong Qiu, or Lee-Ping Wang.
We are also planning to rewrite the ForceBalance documentation with updated instructions and more detailed explanations.
