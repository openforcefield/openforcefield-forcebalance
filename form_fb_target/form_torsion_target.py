# download torsiondrive data from QCArchive server and format as a ForceBalance target
# Author: qyd

import os
import numpy as np
import qcfractal.interface as portal
from forcebalance.molecule import Molecule, bohr2ang

def write_torsion_target(procedure, target_name, removing_unreal=False):
    """
    Fetch data from a QCArchive procedure and write as a ForceBalance Target folder

    Input:
    ------
    procedure: qcfractal.interface.orm.torsiondrive_orm.TorsionDriveORM
        An object that represents a computation

    target_name: string
        Folder name to write the ForceBalance target

    removing_unreal: bool, default False
        If true, remove the atoms that are not "real"
    """
    # check if the target folder name exists
    if os.path.exists(target_name):
        raise RuntimeError(f"Target folder {target_name} already exists")
    # form a Molecule object from the first torsion grid data
    grid_molecules = procedure.final_molecules()
    grid_energies = procedure.final_energies()
    data = next(iter(grid_molecules.values()))
    elem = data['symbols']
    bonds = [(b[0], b[1]) for b in data['connectivity']]
    # remove the "unreal" atoms?
    removing_unreal = removing_unreal and (not all(data['real']))
    if removing_unreal:
        elem = [e for e,f in zip(elem, data['real']) if f is True]
        bonds = [(b[0], b[1]) for b in bonds if data['real'][b[0]] and data['real'][b[1]]]
    # form the Molecule object
    m = Molecule()
    m.Data = {'elem': elem, 'bonds': bonds, 'name': data['name'], 'xyzs': [], 'qm_energies': [], 'comms': []}
    # load all frames
    for grid_id, data in grid_molecules.items():
        geo = np.array(data['geometry']).reshape(-1, 3) * bohr2ang
        if removing_unreal:
            geo = geo[data['real']]
        m.Data['xyzs'].append(geo)
        energy = grid_energies[grid_id]
        m.Data['qm_energies'].append(energy)
        m.Data['comms'].append(f'{data["name"]} at torsion grid {grid_id}')
    # write the data
    os.mkdir(target_name)
    os.chdir(target_name)
    m.write('qdata.txt')
    m.write('geo.xyz')
    os.chdir('..')

def main():
    # connect to server and get result of a torsiondrive run
    client = portal.FractalClient("localhost:7777", verify=False)
    td = client.get_procedures({"procedure": "torsiondrive"})
    # pick the first torsiondrive run from the database for testing
    procedure = td[0]
    # call the above function to make the target
    write_torsion_target(procedure, 'torsion_test')

if __name__ == '__main__':
    main()
