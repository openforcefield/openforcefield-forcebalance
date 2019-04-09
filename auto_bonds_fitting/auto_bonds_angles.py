#!/usr/bin/env python

import os
import copy
import time
import collections
import yaml
import numpy as np

from forcebalance.molecule import Molecule, Elements, bohr2ang
import qcportal as ptl


class FBTargetBuilder:


    def __init__(self, mol2_file, client_conf_file, scan_conf_file=None):
        self.m = Molecule(mol2_file)
        self.qc_mol = self.fb_molecule_to_qc_molecule(self.m)
        # create a client from config file
        self.client = ptl.FractalClient.from_file(client_conf_file)
        # load scan config from scan_conf_file
        self.scan_conf = self.load_scan_conf(scan_conf_file)
        # create output folder
        self.out_folder = os.path.realpath('targets')
        if os.path.exists('targets'):
            raise OSError("Folder targets/ already exist. Please delete the prevous one")
        os.mkdir(self.out_folder)

    def load_scan_conf(self, filename=None):
        """
        Get the scan configuration from a yaml file

        Parameters
        ----------
        filename: str or None
        The input scan config filename (in yaml format). If None, default conf will be used.

        Returns
        -------
        scan_conf: dict
        """
        default_scan_conf = {
            'qm_method': 'HF',
            'qm_basis': 'sto-3g',
            'bond_steps': [-0.2, -0.1, 0.0, 0.1, 0.2],
            'angle_steps': [-20, -10, 0, 10, 20],
        }
        if filename is None:
            return copy.deepcopy(default_scan_conf)
        with open(filename) as infile:
            conf = yaml.load(infile)
        # convert keys to lower case
        conf = {k.lower():v for k,v in conf.items()}
        # check redundant and missing keys
        diff1 = default_scan_conf.keys() - conf.keys()
        if diff1:
            raise ValueError(f"Keys missing in scan_config file {filename}:\n {diff1}")
        diff2 = conf.keys() - default_scan_conf.keys()
        if diff2:
            print(f"Warning: Keys in scan_config file {filename} are ignored:\n {diff2}")
        return conf

    def fb_molecule_to_qc_molecule(self, fb_molecule):
        """ Convert an forcebalance.molecule.Molecule object to a qcportal.Molecule object"""
        e_idxs = [Elements.index(i) for i in fb_molecule.elem]
        coords = fb_molecule.xyzs[0]
        moldata = [[ei] + coord.tolist() for ei, coord in zip(e_idxs, coords)]
        return ptl.Molecule.from_data(moldata, dtype="numpy", units="angstrom")

    def qc_molecule_to_fb_molecule(self, qc_molecule):
        """ Convert an qcportal.Molecule object to a forcebalance.molecule.Molecule object"""
        m = Molecule()
        m.elem = [Elements[i] for i in qc_molecule.atomic_numbers]
        m.xyzs = [qc_molecule.geometry * bohr2ang]
        return m

    def master(self):

        # submit initial optimization
        job_id = self.submit_single_optimization(self.qc_mol)
        self.wait_jobs([job_id])
        self.qc_mol = self.get_optimized_molecule(job_id)

        # submit bond streching jobs
        grid_opt_jobs = self.submit_bond_streching_jobs()

        # submit angle bending jobs
        grid_opt_jobs += self.submit_angle_bending_jobs()

        # submit vibrational hessian jobs
        hessian_job = self.submit_vib_hessian_jobs()

        # wait for grid_opt_jobs to finish
        self.wait_jobs(grid_opt_jobs)

        # collect bond streching job data
        gresult = self.get_grid_optimiztion_results(grid_opt_jobs)
        # flatten the result dict into a list of {'energy': xxx, 'molecule': xxx} records
        flat_records = [record for job_res in gresult.values() for record in job_res.values()]
        self.write_fb_target_abinitio(flat_records)


        # wait for hessian job to finish
        self.wait_jobs([hessian_job], jobtype='compute')
        # collect hessian job data
        hessian_result = self.get_hessian_result(hessian_job)
        self.write_fb_target_hessian(hessian_result)

        # finish




    def submit_single_optimization(self, qc_mol):
        # submit a single optimization job for a qc molecule
        # opt_schema = {}
        # jobId = self.client.add_procedure('optimization', opt_schema)
        # return jobId
        options = {
            "keywords": None,
            "qc_spec": {
                "driver": "gradient",
                "method": self.scan_conf['qm_method'],
                "basis": self.scan_conf['qm_basis'],
                "keywords": None,
                "program": "psi4"
            },
        }
        print(f"Submitting 1 initial optimization job")
        mol_ret = self.client.add_molecules([self.qc_mol])
        r = self.client.add_procedure("optimization", "geometric", options, mol_ret)
        assert len(r.ids) == 1
        return r.ids[0]

    def submit_bond_streching_jobs(self):
        mol_id = self.client.add_molecules([self.qc_mol])[0]
        grid_opt_option_template = {
            "keywords": {
                "preoptimization": False,
                "scans": [{
                    "type": "distance",
                    "indices": None, # To be filled
                    "steps": None, # To be filled
                    "step_type": "relative"
                }]
            },
            "optimization_spec": {
                "program": "geometric",
                "keywords": {
                    "coordsys": "tric",
                    "enforce": 0.1,
                }
            },
            "qc_spec": {
                "driver": "gradient",
                "method": self.scan_conf['qm_method'],
                "basis": self.scan_conf['qm_basis'],
                "keywords": None,
                "program": "psi4",
            },
            "initial_molecule": mol_id,
        }
        all_job_options = []
        strech_steps = self.scan_conf['bond_steps']
        for bond in self.m.bonds:
            job_option = copy.deepcopy(grid_opt_option_template)
            job_option['keywords']['scans'][0]['indices'] = list(bond)
            job_option['keywords']['scans'][0]['steps'] = self.scan_conf['bond_steps']
            all_job_options.append(job_option)
        print(f"Submitting {len(all_job_options)} bond streching grid opt jobs")
        r = self.client.add_service(all_job_options)
        return r.ids

    def submit_angle_bending_jobs(self):
        mol_id = self.client.add_molecules([self.qc_mol])[0]
        grid_opt_option_template = {
            "keywords": {
                "preoptimization": False,
                "scans": [{
                    "type": "angle",
                    "indices": None, # To be filled
                    "steps": None, # To be filled
                    "step_type": "relative"
                }]
            },
            "optimization_spec": {
                "program": "geometric",
                "keywords": {
                    "coordsys": "tric",
                    "enforce": 0.1,
                }
            },
            "qc_spec": {
                "driver": "gradient",
                "method": self.scan_conf['qm_method'],
                "basis": self.scan_conf['qm_basis'],
                "keywords": None,
                "program": "psi4",
            },
            "initial_molecule": mol_id,
        }
        angles = self.m.find_angles()
        #angle_bend_steps = [v/180*3.14159 for v in angle_bend_steps]
        all_job_options = []
        for angle in angles:
            job_option = copy.deepcopy(grid_opt_option_template)
            job_option['keywords']['scans'][0]['indices'] = list(angle)
            job_option['keywords']['scans'][0]['steps'] = self.scan_conf['angle_steps']
            all_job_options.append(job_option)
        print(f"Submitting {len(all_job_options)} angle bending grid opt jobs")
        r = self.client.add_service(all_job_options)
        return r.ids

    def submit_vib_hessian_jobs(self):
        mol_id = self.client.add_molecules([self.qc_mol])[0]
        # submit a hessian job to the server
        method = self.scan_conf['qm_method']
        basis = self.scan_conf['qm_basis']
        r = self.client.add_compute("psi4", method, basis, "hessian", None, mol_id)
        assert len(r.ids) == 1
        return r.ids[0]


    def get_optimized_molecule(self, job_id):
        # get the optimized molecule from a finished job
        qr = self.client.query_procedures(id=job_id)[0]
        return qr.get_final_molecule()

    def get_grid_optimiztion_results(self, grid_opt_jobs):
        """
        Get results of a list of grid optimization jobs
        Return a dictionary in this format:
        {
            grid_opt_job_1: {
                grid_id_1: {
                    'energy': -140.41241,
                    'molecule': qcMol1,
                },
                grid_id_2: {
                    'energy': -140.23241,
                    'molecule': qcMol2,
                }
            },
            grid_opt_job_2: { ... },
            ...
        }
        """
        res = {}
        for job_id in grid_opt_jobs:
            qr = self.client.query_procedures(id=job_id)[0]
            assert qr.status == 'COMPLETE', f'Job {job_id} should be complete, but it is {qr.status.value}'
            # get final energies and geometries for each grid
            energy_dict = qr.get_final_energies()
            molecule_dict = qr.get_final_molecules()
            assert set(energy_dict) == set(molecule_dict), "Keys of energy_dict and molecule_dict should be the same"
            res[job_id] = {}
            for key in energy_dict:
                res[job_id][key] = {
                    'energy': energy_dict[key],
                    'molecule': molecule_dict[key],
                }
        return res

    def get_hessian_result(self, hessian_job_id):
        """
        Get the data from a hessian job
        """
        qr = self.client.query_results(id=hessian_job_id)[0]
        hessian = np.array(qr.return_result, dtype=float)
        # reshape hessian into a matrix, also check the dimensions
        n_of_a = self.m.na
        hessian = hessian.reshape(n_of_a*3, n_of_a*3)
        # get the qcMol
        qcmol = self.client.query_molecules(id=qr.molecule)[0]
        # return a dictionary
        return {'molecule': qcmol, 'hessian': hessian}


    def write_fb_target_abinitio(self, records):
        """ Write a list of {'energy': xxx, 'molecule': xxx, 'name': xxx} records into a new target folder """
        # prepare folder for writing
        target_name = 'abinitio_bond_angles'
        target_folder = os.path.join(self.out_folder, target_name)
        os.mkdir(target_folder)
        os.chdir(target_folder)
        # load data into a fb Molecule
        out_m = Molecule()
        out_m.elem = self.m.elem.copy()
        out_m.xyzs = []
        out_m.qm_energies = []
        out_m.comms = []
        for record in records:
            qcmol = record['molecule']
            energy = record['energy']
            name = record.get('name', 'created by FBTargetBuilder')
            m = self.qc_molecule_to_fb_molecule(qcmol)
            assert m.elem == out_m.elem, 'Elements list of resulting qcmol is not consistent with self.m'
            # append geometry
            out_m.xyzs.append(m.xyzs[0])
            # append energy
            out_m.qm_energies.append(energy)
            # append name
            out_m.comms.append(name)
        # write output
        print(f"Writing {len(records)} frames into targets/abinitio_bond_angles/traj.xyz")
        out_m.write('traj.xyz')
        print(f"Writing {len(records)} frames into targets/abinitio_bond_angles/qdata.txt")
        out_m.write('qdata.txt')

    def write_fb_target_hessian(self, record):
        # prepare folder for writing
        target_name = 'abinitio_hessian'
        target_folder = os.path.join(self.out_folder, target_name)
        os.mkdir(target_folder)
        os.chdir(target_folder)
        # load data into a fb Molecule
        qcmol = record['molecule']
        out_m = self.qc_molecule_to_fb_molecule(qcmol)
        out_m.write('geo.xyz')
        # write hessian matrix
        hessian = record['hessian']
        np.save('hessian', hessian)

    def wait_jobs(self, job_ids, jobtype='procedure', time_interval=5, verbose=True):
        assert jobtype in ['compute', 'procedure']
        while True:
            d_status = collections.defaultdict(int)
            for job_id in job_ids:
                if jobtype == 'procedure':
                    r = self.client.query_procedures(id=job_id)[0]
                elif jobtype == 'compute':
                    r = self.client.query_results(id=job_id)[0]
                status = r.status.value # get string value from RecordStatusEnum
                if r.status == 'ERROR':
                    print(f"Error found in job {jid}")
                    err = r.get_error()
                    if err is not None:
                        print("Error message:")
                        print(err.err_message)
                d_status[status] += 1
            if verbose:
                print(' | '.join(f'{status}:{d_status[status]}' for status in d_status))
            # check if all jobs finished
            if d_status['COMPLETE'] == len(job_ids):
                break
            else:
                time.sleep(time_interval)
        print()

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help='Input mol2 file for a single molecule')
    parser.add_argument("-s", "--scan_config", help='File containing configuration of QM scans')
    parser.add_argument("-c", "--client_config", default='qcportal_config.yaml', help='File containing configuration of QCFractal Client')
    args = parser.parse_args()

    builder = FBTargetBuilder(args.infile, args.client_config, scan_conf_file=args.scan_config)
    builder.master()

if __name__ == '__main__':
    main()