import os, json

def read_optgeo_notes(targets_folder):
    dic = {}
    for fol in os.listdir(targets_folder):
        p = os.path.join(targets_folder, fol)
        if not os.path.isdir(p): continue
        if fol.startswith('optgeo_'):
            #  parse all smiles from each optgeo target
            with open(os.path.join(p,'notes.txt')) as f:
                for line in f:
                    ls = line.split()
                    if len(ls) > 0 and ls[1] == ':':
                        mol2_fnm = os.path.join(p, ls[0]+'.mol2')
                        short_mol2_fnm = 'targets/' + mol2_fnm.split('/')[-2] + '/' + mol2_fnm.split('/')[-1]
                        label = ls[2]
                        smiles = label.rsplit('-',maxsplit=1)[0]
                        mol_id = ls[4]
                        dic[short_mol2_fnm] = {
                            'label' : label,
                            'smiles' : smiles,
                            'mol_id' : mol_id}
    return dic

def main():
    import argparse
    import pickle
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--targets_folder', default='targets')
    args = parser.parse_args()

    dic = read_optgeo_notes(args.targets_folder)

    with open('optgeo_target_mol_info.p', 'wb') as pfile:
        pickle.dump(dic, pfile)

if __name__=='__main__':
    main()
