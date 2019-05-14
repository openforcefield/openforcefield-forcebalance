# forcebalance-qcarchive
Interface between ForceBalance and QCArchive

## Setup guide for a QCArchive local database

References: 
- https://qcfractal.readthedocs.io/en/latest/setup_server.html
- https://github.com/MolSSI/QCFractal/tree/master/examples/local_dataset

1. Install anaconda
    - https://www.anaconda.com/download/#linux

2. Create a new conda environment and activate it
    ```
    conda create -n qcf
    conda activate qcf
    ```

3. Install `qcfractal` and related packages for server
    ```
    conda install -c conda-forge qcfractal qcportal qcengine
    conda install -c psi4 psi4>1.3.1 dftd3
    ```

4. Install ForceBalance
    ```
    conda install -c omnia forcebalance
    ```

5. Start a MongoDB instance
    ```
    MONGOPATH=/tmp/example
    mkdir -p $MONGOPATH
    mongod --dbpath $MONGOPATH
    ```

5. Launch the QCFractal server with a local manager:
    ```
    qcfractal-server mydb --local-manager
    ```

To this point, a server is boot up and jobs can be submitted to it.

Additional steps to build valence or torsiondrive datasets can be found in subfolders.
