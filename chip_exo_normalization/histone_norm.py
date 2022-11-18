import pandas as pd
import argparse

def normalize_cdt(cdt_in, scaling_cdt):
    return (cdt_in.T / scaling_cdt.loc[cdt_in.index].sum(axis=1)).T

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input_cdt')
    parser.add_argument('scaling_cdt')
    parser.add_argument('--out', '-o', default='histone_normalized.cdt', dest='out')
    args = parser.parse_args()

    cdt_in = pd.read_csv(args.input_cdt, sep='\t', index_col=(0, 1))
    scaling_cdt = pd.read_csv(args.scaling_cdt, sep='\t', index_col=(0, 1))
    normalize_cdt(cdt_in, scaling_cdt).to_csv(args.out, sep='\t')
