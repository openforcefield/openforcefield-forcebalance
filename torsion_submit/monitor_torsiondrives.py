#!/usr/bin/env python

import os
import time
import json
from collections import Counter
import qcfractal.interface as ptl

class TorsionMonitor:
    def __init__(self, checkpoint_file, client_conf_file=None, out_folder='td_results'):
        self.load_checkpoint(checkpoint_file)
        self.client = ptl.FractalClient.from_file(client_conf_file) if client_conf_file is not None else None
        self.out_folder = os.path.realpath(out_folder)


    def load_checkpoint(self, filename, verbose=True):
        with open(filename) as infile:
            state = json.load(infile)
        td_jobs = []
        for fname in state:
            if fname == 'scan_conf': continue
            mol_name = os.path.splitext(os.path.basename(fname))[0]
            dihedrals = state[fname]['dihedrals']
            for d in dihedrals:
                job = {
                    'name': mol_name + '_' + d,
                    'mol_name': mol_name,
                    'status': dihedrals[d]['status'],
                    'id': dihedrals[d].get('jobid'),
                }
                td_jobs.append(job)
        self.td_jobs = td_jobs
        if verbose:
            print(f'{len(td_jobs)} jobs found in {filename}')
            self.print_status()

    def get_update(self):
        d_id_jobs = {job['id']: job for job in self.td_jobs if job['id']}
        print(f"Getting updates for {len(d_id_jobs)} jobs by their ids")
        for record in self.client.query_procedures(id=list(d_id_jobs)):
            d_id_jobs[record.id]['status'] = record.status.value

    def print_status(self):
        print('Current Status')
        for status, n in Counter([j['status'] for j in self.td_jobs]).items():
            print(f"-- {status:15s}: {n:5d}")

    def download_complete(self):
        if not os.path.exists(self.out_folder):
            os.mkdir(self.out_folder)
        d_complete_id_jobs = {job['id']:job for job in self.td_jobs if job['status'] == 'COMPLETE'}
        print(f"Downloading results of {len(d_complete_id_jobs)} complete jobs")
        for record in self.client.query_procedures(id=list(d_complete_id_jobs)):
            job = d_complete_id_jobs[record.id]
            folder = os.path.join(self.out_folder, job['mol_name'])
            if not os.path.exists(folder): os.mkdir(folder)
            xyz_filename = os.path.join(folder, job['name'] + '.xyz')
            final_energy_dict = record.get_final_energies()
            final_molecules = record.get_final_molecules()
            with open(xyz_filename, 'w') as xyzfile:
                for grid_id, grid_mol in final_molecules.items():
                    energy = final_energy_dict[grid_id]
                    xyz_str = self.get_xyz_str(grid_mol, title = f"{job['name']} {grid_id} energy = {energy:15.7f}")
                    xyzfile.write(xyz_str + '\n')


    def get_xyz_str(self, qc_mol, title=''):
        elem_list = qc_mol.symbols
        geo = qc_mol.geometry
        noa = len(elem_list)
        lines = [f'{noa}',f'{title}']
        for e, (x,y,z) in zip(elem_list, geo):
            lines.append(f'{e:7s} {x:13.7f} {y:13.7f} {z:13.7f}')
        return '\n'.join(lines)






def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--checkpoint", default="torsion_submit_checkpoint.json", help='Checkpoint file for previous submissions')
    parser.add_argument("-c", "--client_config", help='File containing configuration of QCFractal Client')
    args = parser.parse_args()

    monitor = TorsionMonitor(checkpoint_file=args.checkpoint, client_conf_file=args.client_config)

    monitor.get_update()

    monitor.print_status()

    monitor.download_complete()

if __name__ == "__main__":
    main()