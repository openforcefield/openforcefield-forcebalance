# openforcefield-forcebalance release-1 branch 

This repository contains scripts, whose details are available in the Parsley 1.0.0 paper.  

###  Directory tree
```
|-- README.md
|-- utils/ 
|-- 1_valence_parameter_fitting/
|   |-- 1_dataset_generation/
|   |   |-- roche_set/
|   |   |   |-- 2019-05-16-Roche-Optimization_Set/
|   |   |   |-- 2019-07-09-OpenFF-Optimization-Set/ 
|   |   |   |-- 2019-05-01-OpenFF-Group1-Torsions/ 
|   |   |-- coverage_set/
|   |   |   |-- 2019-06-25-smirnoff99Frost-coverage/
|   |   |   |-- 2019-07-01-smirnoff99Frost-coverage-torsion/ 
|   |-- 2_forcebalance_optimization/
|   |   |-- 1_input_components/ 
|   |   |   |-- optimize.in
|   |   |   |-- forcefield/
|-- 2_benchmarking/
|   |-- test_dataset_generation/ 
|   |   |-- divide_sets.ipynb
|   |   |-- 2019-07-05-eMolecules-force-field-discrepancies-1/
|   |   |-- 2019-09-07-Pfizer-discrepancy-optimization-dataset-1/
|   |   |-- 2019-09-08-fda-optimization-dataset-1/
|-- 3_physical_property_datasets/
|   |-- physical_properties/
```
### Descriptions 

- `utils/README.md`: OpenFF ForceBalance fitting tutorial

- `2019-05-16-Roche-Optimization_Set/`: Roche optimization dataset construction and submission [(original location)](https://github.com/openforcefield/qca-dataset-submission/tree/master/2019-05-16-Roche-Optimization_Set)

- `2019-07-09-OpenFF-Optimization-Set/`: Roche hessian dataset construction and submission [(original location)](https://github.com/openforcefield/qca-dataset-submission/tree/master/2019-07-09-OpenFF-Optimization-Set)

- `2019-05-01-OpenFF-Group1-Torsions/`: Roche torsiondrive dataset construction and submission [(original location)](https://github.com/openforcefield/qca-dataset-submission/tree/master/2019-05-01-OpenFF-Group1-Torsions)

- `019-06-25-smirnoff99Frost-coverage/`: coverage optimization and hessian dataset construction and submission [(original location)](https://github.com/openforcefield/qca-dataset-submission/tree/master/2019-06-25-smirnoff99Frost-coverage)

- `2019-07-01-smirnoff99Frost-coverage-torsion/`: coverage torsiondrive dataset construction and submission [(original location)](https://github.com/openforcefield/qca-dataset-submission/tree/master/2019-07-01-smirnoff99Frost-coverage-torsion)

- `optimize.in`: input file for ForceBalance optimization [(original location)](https://github.com/openforcefield/openforcefield-forcebalance/releases/tag/v1.0.0-RC2)

- `forcefield/`: input forcefield directory for ForceBalance optimizaiton [(original location)](https://github.com/openforcefield/openforcefield-forcebalance/releases/tag/v1.0.0-RC2)

- `divide_sets.ipynb`: script for dividing a large benchmark set into these two sets using path-based fingerprint method [(original location)](https://github.com/openforcefield/release-1-benchmarking/blob/master/QM_molecule_selection/divide_sets.ipynb)

- `2019-07-05-eMolecules-force-field-discrepancies-1/`: eMoleculesa discrepancy set optimization dataset construction and submission [(original location)](https://github.com/openforcefield/qca-dataset-submission/blob/master/2019-07-05%20eMolecules%20force%20field%20discrepancies%201)

- `2019-09-07-Pfizer-discrepancy-optimization-dataset-1/`: Pfizer discrepancy set optimization dataset construction and submission [(original location)](https://github.com/openforcefield/qca-dataset-submission/tree/master/2019-09-07-Pfizer-discrepancy-optimization-dataset-1)

- `2019-09-08-fda-optimization-dataset-1/`: FDA optimization dataset construction and submission [(original location)](https://github.com/openforcefield/qca-dataset-submission/tree/master/2019-09-08-fda-optimization-dataset-1)

- `physical_properties/` : contains scripts for automated curation of physical property datasets and benchmarking results against condensed phase liquid properties [(original location)](https://github.com/openforcefield/release-1-benchmarking/tree/master/physical_properties)