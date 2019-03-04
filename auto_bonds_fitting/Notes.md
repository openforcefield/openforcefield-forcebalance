# Notes for automated fitting using QCFractal and ForceBalance


1. Start with mol2 file with topology

2. Create a `qcfractal.interface.Molecule` object for the target molecule
    ```
    hooh = portal.Molecule({"symbols": ["H", "O", "O", "H"], "geometry": [0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4]})
    client.add_molecules({"hooh": hooh})
    ```
    We should optimize the geometry first.

3. (Skipped) For fitting bond and angle, one option is to run an AIMD simulation with high temperature, save the trajectory. (Not available on QCFractal yet) Add each frame of the trajectory as a `Molecule` in the database, then submit them by `add_compute` or `add_procedure` or `add_service`.

4. Another option is to use "grid optimization". The user need to select each bond and perturb them using `add_service` for grid optimization. (Doc/Example missing) Then each angle can be scanned. (It supports coupled scan if needed) (We should also use Hessians)
    ```
    service = GridOptimizationInput(**{
        "gridoptimization_meta": {
            "starting_grid":
            "zero",
            "scans": [{
                "type": "distance",
                "indices": [1, 2],
                "steps": [1.0, 1.1] # 9-11 points in real calculation
            }]
        },
        "optimization_meta": {
            "program": "geometric",
            "coordsys": "tric",
        },
        "qc_meta": {
            "driver": "gradient",
            "method": "UFF",
            "basis": "",
            "options": None,
            "program": "rdkit",
        },
        "initial_molecule": [hooh],
    })

    ret = client.add_service(service)
    ```
    https://github.com/MolSSI/QCFractal/issues/112


5. For fitting torsion parameters, user need to pick the interested torsions and their scan options, then submit a `TorsionDrive` by `add_service`
    ```
    tdinput = {
        "initial_molecule": [hooh],
        "torsiondrive_meta": {
            "dihedrals": [[0, 1, 2, 3]],
            "grid_spacing": [90]
        },
        "optimization_meta": {
            "program": "geometric",
            "coordsys": "tric",
        },
        "qc_meta": {
            "driver": "gradient",
            "method": "UFF",
            "basis": None,
            "options": None,
            "program": "rdkit",
        },
    }

    ret = client.add_service(tdinput)
    ```
    https://github.com/MolSSI/QCFractal/blob/master/examples/parsl_torsiondrive/compute_torsion.py


6. For fitting non-bonded parameters, possible to use interaction energy `Dataset`:
    https://qcfractal.readthedocs.io/en/latest/collection-dataset.html


7. Collected the results of each above job, and save them on disk as ForceBalance targets. (Tools ready for torsiondrive results, should be straight forward for grid optimization)

8. Prepare ForceBalance input file and force field file, run ForceBalance to get fitted parameters.

9. Evaluate fit quality, by checking QM vs MM energy and parameter changes. Prepare new fitting targets and repeat above steps if needed.