#!/usr/bin/env python

import os
import sys
import copy
import yaml
import json
import numpy as np

from forcebalance.molecule import Molecule, Elements, bohr2ang
import qcfractal.interface as ptl

from find_dihedrals import DihedralSelector

class TorsionSubmitter:
    def __init__(self, scan_conf_file=None, client_conf_file=None):
        # load scan config from scan_conf_file
        self.scan_conf = self.load_scan_conf(scan_conf_file)
        # create a client from config file
        self.client = ptl.FractalClient.from_file(client_conf_file) if client_conf_file is not None else None
        # load checkpoint file
        self.checkpoint_filename = "torsion_submit_checkpoint.json"
        self.load_checkpoint(self.checkpoint_filename)
        # sync scan config with checkpoint
        self.sync_scan_conf()
        # hold all json option to submit
        self.json_submit_options = []

    def load_checkpoint(self, filename):
        if os.path.isfile(filename):
            with open(filename) as infile:
                self.state = json.load(infile)
        else:
            print(f"checkpoint file {filename} not found, starting with empty state")
            self.state = {}

    def write_checkpoint(self):
        with open(self.checkpoint_filename, 'w') as outfile:
            json.dump(self.state, outfile, indent=2)

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
            'qm_method': 'b3lyp',
            'qm_basis': '6-31G*',
            'grid_spacing': 15,
            'energy_upper_limit': 0.05,
        }
        if filename is None:
            print("scan_configure file not provided, using default config")
            fn = "td_scan_configure.yaml"
            conf = copy.deepcopy(default_scan_conf)
            with open(fn, 'w') as outfile:
                yaml.dump(conf, outfile, default_flow_style=False)
            print(f"scan configure saved as {fn}")
            return conf
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

    def sync_scan_conf(self):
        if 'scan_conf' in self.state:
            existing_conf = self.state['scan_conf']
            current_conf = self.scan_conf
            # make sure the current configuration matches the previous one
            assert json.dumps(existing_conf, sort_keys=True) == json.dumps(current_conf, sort_keys=True), \
                f'Error: previous scan conf {existing_conf} not consistent with current conf {current_conf}'
        else:
            self.state['scan_conf'] = self.scan_conf

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

    def submit_molecule(self, filename, to_json=False):
        print(f"\n*** Submitting torsion scans for {filename} ***")
        m = Molecule(filename)
        qc_mol = self.fb_molecule_to_qc_molecule(m)
        if to_json is False:
            mol_id = self.client.add_molecules([qc_mol])[0]
        else:
            mol_id = qc_mol.json_dict()
        dihedral_selector = DihedralSelector(filename)
        dihedral_list = dihedral_selector.find_dihedrals(dihedral_filter='heavy_no_ring')
        all_job_options = []
        for dihedral in dihedral_list:
            torsiondrive_options = {
                "initial_molecule": mol_id,
                "keywords": {
                    "dihedrals": [dihedral],
                    "grid_spacing": [self.scan_conf['grid_spacing']],
                    "energy_upper_limit": self.scan_conf['energy_upper_limit'],
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
            }
            all_job_options.append(torsiondrive_options)
        if to_json is False:
            print(f"Submitting {len(all_job_options)} torsiondrive jobs")
            r = self.client.add_service(all_job_options)
            # store the state in checkpoint
            self.state[filename] = {}
            for dihedral, jobid in zip(dihedral_list, r.ids):
                self.state[filename]['-'.join(dihedral)] = {'jobid': jobid, 'status': 'submitted'}
        else:
            print(f"Saving {len(all_job_options)} torsiondrive options to json")
            self.json_submit_options.extend(all_job_options)

    def write_submitted_json(self, filename):
        with open(filename, 'w') as outfile:
            json.dump(self.json_submit_options, outfile, indent=2)
            print(f'Total {len(self.json_submit_options)} torsiondrive options written to {filename}')

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("infiles", nargs='*', help='Input mol2 file for a single molecule')
    parser.add_argument("-s", "--scan_config", help='File containing configuration of QM scans')
    parser.add_argument("-c", "--client_config", help='File containing configuration of QCFractal Client')
    parser.add_argument("-j", "--save_json", action="store_true", default=False, help='If specified, will not submit but save all submit job in a json.')
    args = parser.parse_args()

    print(' '.join(sys.argv))

    submitter = TorsionSubmitter(scan_conf_file=args.scan_config, client_conf_file=args.client_config)
    for f in args.infiles:
        submitter.submit_molecule(f, to_json=args.save_json)

    if args.save_json:
        submitter.write_submitted_json("submit_torsion_options.json")
    else:
        submitter.write_checkpoint()


if __name__ == '__main__':
    main()
