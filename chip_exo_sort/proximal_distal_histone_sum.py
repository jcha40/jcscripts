import pandas as pd
import argparse

def get_sort(sense_cdts, anti_cdts, proximal_idx=slice(400, 551), distal_idx=slice(449, 600)):
    idx = sense_cdts[0].index
    sort_ser = pd.Series(dict.fromkeys(idx, 0))
    for sense_cdt, anti_cdt in zip(sense_cdts, anti_cdts):
        sort_ser += sense_cdt.iloc[:, proximal_idx].sum(axis='columns') + anti_cdt.iloc[:, distal_idx].sum(axis='columns')
    return sort_ser.sort_values(ascending=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sense-cdt', '-s', nargs='+', dest='sense_cdt')
    parser.add_argument('--anti-cdt', '-a', nargs='+', dest='anti_cdt')
    parser.add_argument('--proximal-idx', '-p', type=int, nargs=2, dest='proximal_idx', default=(400, 551))
    parser.add_argument('--distal-idx', '-d', type=int, nargs=2, dest='distal_idx', default=(449, 600))
    parser.add_argument('--out', '-o', default='proximal_distal_histone_sum_sort.tsv', dest='out')
    args = parser.parse_args()

    sense_cdt = pd.read_csv(args.sense_cdt, sep='\t', index_col=(0, 1))
    anti_cdt = pd.read_csv(args.anti_cdt, sep='\t', index_col=(0, 1))
    proximal_idx = slice(*args.proximal_idx)
    distal_idx = slice(*args.distal_idx)
    get_sort(sense_cdt, anti_cdt, proximal_idx, distal_idx).to_csv(args.out, sep='\t', header=None)
