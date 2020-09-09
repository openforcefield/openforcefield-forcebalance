import os
import numpy as np

def read_abinitio_EnergyCompare(tmp_folder):
    dic = {}
    for fol in os.listdir(tmp_folder):
        p = os.path.join(tmp_folder, fol)
        if not os.path.isdir(p): continue
        if fol.startswith('abinitio_'):
            Efile = os.path.join(p,'iter_0000/EnergyCompare.txt')
            data = np.loadtxt(Efile, ndmin=2)
            qme = data[:,0]
            mme = data[:,1]
            delta = data[:,2]
            avg_delta = np.mean(abs(np.array(delta)))
            dic[fol] = avg_delta
    return dic

def read_abinitio_notes(targets_folder):
    dic = {}
    for fol in os.listdir(targets_folder):
        p = os.path.join(targets_folder, fol)
        if not os.path.isdir(p): continue
        if fol.startswith('abinitio_'):
            #  parse all smiles from each optgeo target
            with open(os.path.join(p,'notes.txt')) as f:
                mol_ids = []
                for line in f:
                    ls = line.split()
                    if len(ls) > 0 and ls[1] == 'molecule_id':
                        label = ls[0]
                        smiles = label.rsplit('-',maxsplit=1)[0]
                        mol_id = ls[2]
                        mol_ids.append(mol_id)
            dic[fol] = {
                'smiles' : smiles,
                'mol_ids' : mol_ids}
    return dic

def sort_print_obj_table(molecule_obj_dict, mol_info):
    sorted_molecules = sorted(molecule_obj_dict, key=lambda x: molecule_obj_dict[x])
    print("Results")
    print(f"idx       {'folder':^70s} {'SMILES':^50s} objective")
    print('-'*120)
    for i, folder in enumerate(sorted_molecules):
        foldernm = 'targets/' + folder
        obj = molecule_obj_dict[folder]
        smiles = mol_info[folder]['smiles']
        print(f'{i:<7}   {foldernm:70s} {smiles:50s} {obj:9.5f}')


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--targets_folder', default='targets')
    parser.add_argument('-tmp', '--tmp_folder')
    args = parser.parse_args()

    mol_info = read_abinitio_notes(args.targets_folder)
    molecule_obj_dict = read_abinitio_EnergyCompare(args.tmp_folder)
    sort_print_obj_table(molecule_obj_dict, mol_info)

if __name__=='__main__':
    main()
