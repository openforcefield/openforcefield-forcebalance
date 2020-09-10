import os
import json

def find_unique_smiles(targets_folder):
    unique_smiles = {}
    dataset_entries = {}
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
                        dataset_entries[label] = p
                        sm = label.rsplit('-',maxsplit=1)[0]
                        unique_smiles[sm] = p
        elif fol.startswith('vibfreq_'):
            # parse all smiles from each vibfreq target
            with open(os.path.join(p, 'note.txt')) as f:
                for line in f:
                    label = line.rsplit(maxsplit=1)[-1]
                    dataset_entries[label] = p
                    sm = label.rsplit('-',maxsplit=1)[0]
                    unique_smiles[sm] = p
        elif fol.startswith('td_'):
            # parse all smiles from each torsionprofile target
            with open(os.path.join(p, 'metadata.json')) as f:
                d = json.load(f)
                label = d["entry_label"]
                dataset_entries[label] = p
                sm = d['attributes']['canonical_isomeric_smiles']
                unique_smiles[sm] = p
    return unique_smiles, dataset_entries


def print_each_line(unique_smiles, **kwargs):
    for s in sorted(unique_smiles):
        print(f"{s:80} {unique_smiles[s]}", **kwargs)

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--targets_folder', default='targets')
    args = parser.parse_args()

    unique_smiles, dataset_entries = find_unique_smiles(args.targets_folder)

    print_each_line(unique_smiles)

    with open("dataset_entries.txt", "w") as f:
        print_each_line(dataset_entries, file=f)

if __name__ == "__main__":
    main()