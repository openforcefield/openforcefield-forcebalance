import os

def find_unique_smiles(targets_folder):
    unique_smiles = {}
    for fol in os.listdir(targets_folder):
        p = os.path.join(targets_folder, fol)
        if not os.path.isdir(p): continue
        if fol.startswith('optgeo_'):
            # parse all smiles from each optgeo target
            with open(os.path.join(p, 'notes.txt')) as f:
                for line in f:
                    ls = line.split()
                    if len(ls) > 0 and ls[1] == ':':
                        label = ls[2]
                        sm = label.rsplit('-',maxsplit=1)[0]
                        unique_smiles[sm] = p
        elif fol.startswith('vibfreq_'):
            # parse all smiles from each vibfreq target
            with open(os.path.join(p, 'note.txt')) as f:
                for line in f:
                    label = line.rsplit(maxsplit=1)[-1]
                    sm = label.rsplit('-',maxsplit=1)[0]
                    unique_smiles[sm] = p
        elif fol.startswith('td_'):
            # parse all smiles from each torsionprofile target
            with open(os.path.join(p, 'note.txt')) as f:
                for line in f:
                    label = line.rsplit(maxsplit=1)[-1]
                    sm = label.replace(':1','').replace(':2','').replace(':3','').replace(':4','')
                    unique_smiles[sm] = p
    return unique_smiles


def print_each_line(unique_smiles):
    for s in sorted(unique_smiles):
        print(f"{s:80} {unique_smiles[s]}")

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--targets_folder', default='targets')
    args = parser.parse_args()

    unique_smiles = find_unique_smiles(args.targets_folder)

    print_each_line(unique_smiles)

if __name__ == "__main__":
    main()