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
        self.downloaded_json_fn = os.path.join(self.out_folder, 'downloaded.json')

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
        print(f"Updating status for {len(d_id_jobs)} jobs by their ids")
        # check the out_folder for jobs that are downloaded already
        if os.path.exists(self.downloaded_json_fn):
            downloaded_jobs = json.load(open(self.downloaded_json_fn))
            for job in downloaded_jobs:
                # update status to "DOWNLOADED"
                d_id_jobs[job['id']]['status'] = "DOWNLOADED"
        # pull status from server for other jobs
        query_job_ids = [job['id'] for job in d_id_jobs.values() if job['status'] != 'DOWNLOADED']
        for record in self.client.query_procedures(id=query_job_ids):
            d_id_jobs[record.id]['status'] = record.status.value

    def print_status(self):
        print('< Current Status >')
        print('|', end='')
        for status, n in Counter([j['status'] for j in self.td_jobs]).items():
            print(f"{status:>15s}: {n:<5d}", end='|')
        print()

    def download_complete(self):
        if not os.path.exists(self.out_folder):
            os.mkdir(self.out_folder)
        d_complete_id_jobs = {job['id']:job for job in self.td_jobs if job['status'] == 'COMPLETE'}
        n = len(d_complete_id_jobs)
        print(f"Downloading results of {n} complete jobs")
        for i, record in enumerate(self.client.query_procedures(id=list(d_complete_id_jobs)), 1):
            job = d_complete_id_jobs[record.id]
            print(f"{i:>3d}/{n:<3d} Downloading results for job {job['name']}")
            # prepare folder
            folder = os.path.join(self.out_folder, job['mol_name'])
            if not os.path.exists(folder): os.mkdir(folder)
            # get data
            final_energy_dict = record.get_final_energies()
            final_molecules = record.get_final_molecules()
            # write xyz file
            xyz_filename = os.path.join(folder, job['name'] + '.xyz')
            with open(xyz_filename, 'w') as xyzfile:
                for grid_id in sorted(final_molecules):
                    grid_mol = final_molecules[grid_id]
                    energy = final_energy_dict[grid_id]
                    xyz_str = self.get_xyz_str(grid_mol, title = f"{job['name']} {grid_id} energy = {energy:15.7f}")
                    xyzfile.write(xyz_str + '\n')
            # save energy curve plot as pdf
            plot_filename = os.path.join(folder, job['name'] + '.pdf')
            self.plot_1d_energies(final_energy_dict, plot_filename, title=job['name'])
            # change status of job
            job['status'] = 'DOWNLOADED'
            job['saved_file'] = os.path.relpath(xyz_filename)
        # save downloaded job information in file
        with open(self.downloaded_json_fn, 'w') as jsonfile:
            downloaded_jobs = [job for job in self.td_jobs if job['status'] == 'DOWNLOADED']
            json.dump(downloaded_jobs, jsonfile, indent=2)

    def get_xyz_str(self, qc_mol, title=''):
        elem_list = qc_mol.symbols
        # convert geometry unit Bohr -> Angstrom
        geo = qc_mol.geometry * 0.529177
        noa = len(elem_list)
        lines = [f'{noa}',f'{title}']
        for e, (x,y,z) in zip(elem_list, geo):
            lines.append(f'{e:7s} {x:13.7f} {y:13.7f} {z:13.7f}')
        return '\n'.join(lines)

    def plot_1d_energies(self, energy_dict, filename, title=''):
        if not energy_dict:
            print("Empty energy dict, skip plotting")
            return
        grid_ids = sorted(energy_dict.keys())
        x_dihedrals = [gid[0] for gid in grid_ids]
        y_energies = [energy_dict[gid] for gid in grid_ids]
        # convert to relative energies in kcal/mol
        import numpy as np
        y_energies = np.array(y_energies)
        y_energies = (y_energies - np.min(y_energies)) * 627.509
        # plot
        import matplotlib.pyplot as plt
        plt.style.use('ggplot')
        plt.Figure()
        plt.plot(x_dihedrals, y_energies, '-o')
        plt.xlabel("Dihedral Angle [degrees]")
        plt.ylabel("Relative Energies [kcal/mol]")
        plt.title(title)
        plt.savefig(filename)
        plt.close()




def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--checkpoint", default="torsion_submit_checkpoint.json", help='Checkpoint file for previous submissions')
    parser.add_argument("-c", "--client_config", default='client_config.yaml', help='File containing configuration of QCFractal Client')
    args = parser.parse_args()

    monitor = TorsionMonitor(checkpoint_file=args.checkpoint, client_conf_file=args.client_config)

    monitor.get_update()

    monitor.print_status()

    monitor.download_complete()

if __name__ == "__main__":
    main()