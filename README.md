# forcebalance-qcarchive
Interface between ForceBalance and QCArchive

## Setup guide for a QCArchive local database
reference: https://github.com/MolSSI/QCFractal/tree/master/examples/local_dataset

1. Install MongoDB via Docker
    - The clean and easy way to install MongoDB for testing is via Docker:
        https://docs.docker.com/install/linux/docker-ce/ubuntu/#install-docker-ce

    + Once installed, use command
        ```
        docker pull mongo
        ```
        to pull the official MongoDB image from DockerHub
        Then run
        ```
        docker run mongo -p 27017:27017 -d
        ```
        to run the container in background with shared default port `27017`

2. Install anaconda
    ```
    https://www.anaconda.com/download/#linux
    ```

3. clone the QCFractal git repository:
    ```
    git clone https://github.com/MolSSI/QCFractal.git
    ```

4. create a new conda env for testing
    ```
    cd QCFractal
    python devtools/scripts/conda_env.py -n=qcf -p=3.6 devtools/conda-envs/openff.yaml
    ```

5. Install QCFractal
    ```
    conda activate qcf
    pip install -e .
    ```

6. Install ForceBalance
    ```
    conda install -c omnia forcebalance
    ```

7. Launch the QCFractal server:
    ```
    qcfractal-server qca_parsl_testing
    ```

8. In a new terminal window, launch the `parsl_manager`:
    ```
    conda activate qcf
    cd QCFractal/examples/parsl_torsiondrive
    python parsl_manager.py
    ```

9. In a new terminal window, submit an example torsiondrive job:
    ```
    conda activate qcf
    cd QCFractal/examples/parsl_torsiondrive
    python compute_torsion.py
    ```
    The other terminal windows of `parsl_manager` and `qcfractal-server` should status of the jobs running. The example job should finish fairly quick.

10. After the job finishes, run the script to pull data from server and form a ForceBalance target:
    ```
    cd forcebalance-qcarchive/form_fb_target
    python form_torsion_target.py
    ```
    A new folder `torsion_test/` should be created, that contains a `qdata.txt` file and an `.xyz` file for ForceBalance to read.

