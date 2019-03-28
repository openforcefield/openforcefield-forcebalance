# coding: utf-8
from auto_bonds_angles import FBTargetBuilder
self = builder = FBTargetBuilder('meoh.mol2', 'qcportal_config.yaml')
job_id = self.submit_single_optimization(self.qc_mol)
self.wait_jobs([job_id])
self.qc_mol = self.get_optimized_molecule(job_id)
self.qc_mol.id
grid_opt_jobs = self.submit_bond_streching_jobs()
grid_opt_jobs += self.submit_angle_bending_jobs()
self.wait_jobs(grid_opt_jobs)
