import os
import pickle
import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


def read_pickle_data(fnm):
    """
    Read data from pickled file generated by plot_optgeo_each_smirks.py

    The content of the pickle file:
    data_save = {
        'ffxml': args.ffxml,
        'tmp_folder': tmp_folder,
        'iter_folder': iter_folder,
        'data_qm_v_mm': data_qm_v_mm,
    }

    The main content is in data_qm_v_mm
    data_qm_v_mm[fftype][sid] = [
        {
            'mol2_fnm': mol2_fnm,
            'atom_indices': atom_indices,
            'smirks': smirks,
            'id': sid,
            'qm': qm,
            'mm': mm,
            'mm_iter0': None or mm of iter0,
        },
        ...
    ]
    """
    print(f'reading saved data from {fnm}')
    with open(fnm, 'rb') as pfile:
        data_save = pickle.load(pfile)
        data_qm_v_mm = data_save['data_qm_v_mm']
    return data_qm_v_mm

def compute_rmse(ref_data_list, target_data_list):
    a = np.array(ref_data_list)
    b = np.array(target_data_list)
    rmse = np.sqrt(np.sum((a-b)**2) / len(a))
    return rmse

def aggregate_compare_rmse_data(optgeo_data_list):
    agg_data = {}
    orig_data = optgeo_data_list[0]
    for fftype in orig_data:
        agg_data[fftype] = []
        for sid in sorted(orig_data[fftype], key=lambda s: int(s[1:] if s[-1].isdigit() else s[1:-1])):
            qm_values = [d['qm'] for d in orig_data[fftype][sid]]
            rmse_list = []
            for data in optgeo_data_list:
                mm_values = [d['mm'] for d in data[fftype][sid]]
                rmse = compute_rmse(qm_values, mm_values)
                rmse_list.append(rmse)
            agg_data[fftype].append({
                'sid': sid,
                'rmse_list': rmse_list,
            })
    return agg_data

def plot_compare_rmse_smirks(agg_data, fnm='compare_optgeo_rmse.pdf', legends=None):
    """ Generate bar plots for comparing benchmark results """
    xlabels = {
        'bonds': 'Bond Length RMSE (Angstrom)',
        'angles': 'Bond Angles RMSE (Degrees)',
        'propertorsions': 'Torsion Angles RMSE (Degrees)',
        'impropertorsions': 'Improper Torsion Angles RMSE (Degrees)'
    }
    with PdfPages(fnm) as pdf:
        for fftype, sid_data_list in agg_data.items():
            sid_list = [d['sid'] for d in sid_data_list]
            # transpose so each row is rmse array for one forcefield
            rmse_value_matrix = np.array([d['rmse_list'] for d in sid_data_list]).T
            # make plot
            n = len(sid_list)
            n_bars = len(rmse_value_matrix)
            y_gap = n_bars*0.3 + 0.2
            y_pos = np.arange(n) * y_gap
            plt.figure(figsize=(8.5, n*0.06*n_bars+1.2))
            for i, rmse_values in enumerate(rmse_value_matrix):
                label = legends[i] if legends is not None else None
                if i == 0:
                    plt.barh(y_pos, rmse_values, tick_label=sid_list, height=0.3, align='center', label=label)
                else:
                    y_shift = 0.3 * i
                    plt.barh(y_pos + y_shift, rmse_values, height=0.3, align='center', label=label)
            # adjust the y range, and invert the yaxis
            plt.ylim(y_pos[-1]+y_gap, y_pos[0]-1)
            # adjust the x range
            xmin = 0
            xmax = rmse_value_matrix.max()
            padding = (xmax - xmin) * 0.01
            plt.xlim(xmin, xmax+padding)
            plt.xlabel(xlabels[fftype.lower()])
            if legends:
                plt.legend()
            # save
            plt.title(f'RMSE comparison for {fftype}')
            plt.tight_layout()
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()
    print(f"RMSE compare plots saved as {fnm}")

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('data_pickle_files', nargs='+', help='pickled data generated by plot_optgeo_each_smirks.py')
    args = parser.parse_args()

    file_names = [os.path.splitext(os.path.basename(f))[0] for f in args.data_pickle_files]
    optgeo_data_list = [read_pickle_data(f) for f in args.data_pickle_files]

    agg_data = aggregate_compare_rmse_data(optgeo_data_list)

    plot_compare_rmse_smirks(agg_data, legends=file_names)

if __name__ == "__main__":
    main()