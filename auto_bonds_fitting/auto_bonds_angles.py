#!/usr/bin/env python

import os
import copy
from forcebalance.molecule import Molecule, Elements, bohr2ang
import qcportal as ptl

class FBTargetBuilder:
    def __init__(self, mol2_file, conf_file):
        self.m = Molecule(mol2_file)
        self.qc_mol = self.fb_molecule_to_qc_molecule(self.m)
        self.client = ptl.FractalClient.from_file(conf_file)
        self.out_folder = 'targets'
        if os.path.exists('targets'):
            raise OSError("Folder targets/ already exist. Please delete the prevous one")
        os.mkdir(self.out_folder)

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

        # create the molecule on server
        self.client.add_molecules([self.qc_mol])

        # submit initial optimization
        job = self.submit_single_optimization(self.qc_mol)
        self.wait_jobs([job])
        self.qc_mol = self.get_optimized_molecule(job)

        # add the optimized molecule to server
        self.client.add_molecules([self.qc_mol])

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
        self.write_fb_target_abinitio(gresult)


        # wait for hessian job to finish
        self.wait_jobs([hessian_job])
        # collect hessian job data
        hresult = self.get_hessian_result(hessian_job)
        self.write_fb_target_hessian(hresult)

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
                "method": "HF",
                "basis": "sto-3g",
                "keywords": None,
                "program": "psi4"
            },
        }
        r = self.client.add_procedure("optimization", "geometric", options, qc_mol)
        assert len(r.ids) == 1
        return r.ids[0]

    def submit_bond_streching_jobs(self):
        grid_opt_option_template = {
            "keywords": {
                "preoptimization": True,
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
                }
            },
            "qc_spec": {
                "driver": "gradient",
                "method": "HF",
                "basis": "sto-3g",
                "keywords": None,
                "program": "psi4",
            },
            "initial_molecule": self.qc_mol,
        }
        all_job_options = []
        # perturb the length of each bond by -30%, -20%, -10%, -5%, 5% ...
        strech_steps = [-0.25, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]
        for bond in self.m.bonds:
            job_option = copy.deepcopy(grid_opt_option_template)
            job_option['keywords']['scans'][0]['indices'] = list(bond)
            job_option['keywords']['scans'][0]['steps'] = strech_steps
            all_job_options.append(job_option)
        r = self.client.add_service(all_job_options)
        return r.ids

    def submit_angle_bending_jobs(self):
        grid_opt_option_template = {
            "keywords": {
                "preoptimization": True,
                "scans": [{
                    "type": "angle",
                    "indices": None, # To be filled
                    "steps": None, # To be filled
                    "step_type": "absolute"
                }]
            },
            "optimization_spec": {
                "program": "geometric",
                "keywords": {
                    "coordsys": "tric",
                }
            },
            "qc_spec": {
                "driver": "gradient",
                "method": "HF",
                "basis": "sto-3g",
                "keywords": None,
                "program": "psi4",
            },
            "initial_molecule": self.qc_mol,
        }
        angles = self.m.find_angles()
        # perturb the value of each angle by -30, -20, -10, 5, 10 ...
        angle_bend_steps = [-30, -20, -10, 0, 10, 20, 30]
        all_job_options = []
        for angle in angles:
            job_option = copy.deepcopy(grid_opt_option_template)
            job_option['keywords']['scans'][0]['indices'] = list(angle)
            job_option['keywords']['scans'][0]['steps'] = angle_bend_steps
            all_job_options.append(job_option)
        r = self.client.add_service(all_job_options)
        return r.ids

    def submit_vib_hessian_jobs(self):
        #schema = {}
        # submit a hessian job to the server
        # job_id = self.client.add_procedule('hessian', schema, self.m)
        # return job_id
        return 1

    def get_optimized_molecule(self, job_id):
        # get the optimized molecule from a finished job
        qc_mol = 1
        return qc_mol

    def wait_jobs(self, job_list):
        # waiting for all jobs in job_list to finish
        return None

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help='Input mol2 file for a single molecule')
    parser.add_argument("-c", "--fractal_config", default='qcportal_config.yaml', help='File containing configuration of QCFractal Client')
    args = parser.parse_args()

    builder = FBTargetBuilder(args.infile, args.fractal_config)
    builder.master()

if __name__ == '__main__':
    main()