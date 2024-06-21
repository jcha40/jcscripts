import numpy as np
import os
import re
import h5py
from filelock import SoftFileLock
import argparse

def pileup(scidx_fn: str, bed_fns: list, chrom_sizes: str, out: str, control: bool = False,
           ref_pattern: re.Pattern = re.compile(r'(.+)_(\d+)bp.bed'),
           sample_pattern: re.Pattern = re.compile(r'(\d+)_(.+?)_i5006_BY4741_-_YPD_(.+?)_XO_FilteredBAM')):
    """
    scidx_fn: str
        Path to the scidx file
    bed_fns: list
        List of paths to bed files
    chrom_sizes: str
        Path to the chrom.sizes file
    out: str
        Path to the output h5 file
    control: bool
        If True, the output will be formatted as a control pileup
    ref_pattern: re.Pattern
        Regular expression pattern to extract the reference name from the bed file name
    sample_pattern: re.Pattern
        Regular expression pattern to extract the replicate, target, and condition from the scidx file name
    """
    sizes = {}
    chrom_dict = {}
    with open(chrom_sizes) as f:
        for line in f:
            chrom, size = line.strip().split()
            sizes[chrom] = int(size)
            chrom_dict[chrom] = np.zeros((sizes[chrom], 2), dtype=np.int32)
    
    with open(scidx_fn, 'r') as f:
        contig = None
        f.readline()
        f.readline()
        for line in f:
            line = line.strip().split('\t')
            contig = line[0]
            pos = int(line[1]) - 1
            chrom_dict[contig][pos, 0] = int(line[2])
            chrom_dict[contig][pos, 1] = int(line[3])
    
    sample_id, target, cond = sample_pattern.match(os.path.basename(scidx_fn)).groups()
    lock = SoftFileLock('{}.lock'.format(out))
    with lock:
        with h5py.File(out, 'a') as h5:
            target_cond = '{}-{}'.format(target, cond)
            if h5.get(target_cond) is None:
                h5.create_group(target_cond)
            target_group = h5[target_cond]
            if target_group.get(sample_id) is None:
                target_group.create_group(sample_id)
            
            if control:
                if h5.get('controls') is None:
                    h5.create_group('controls')
                controls = set(h5['controls'][:])
                controls.add(target_cond)
                h5['controls'] = list(controls)

    for bed_fn in bed_fns:
        print(bed_fn)
        ref_name, length = ref_pattern.match(os.path.basename(bed_fn)).groups()
        length = int(length)
        forward_comp = np.zeros(length, dtype=np.int32)
        reverse_comp = np.zeros(length, dtype=np.int32)
        with open(bed_fn, 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                contig = line[0]
                start = int(line[1])
                end = int(line[2])
                if start < 1 or end > sizes[contig] - 1:
                    print('Skipping {}:{}-{} (out of bounds)'.format(contig, start, end))
                    continue
                if line[5] == '-':
                    forward_comp += chrom_dict[contig][end - 1:start - 1:-1, 1]
                    reverse_comp += chrom_dict[contig][end - 1:start - 1:-1, 0]
                else:
                    forward_comp += chrom_dict[contig][start:end, 0]
                    reverse_comp += chrom_dict[contig][start:end, 1]
        
        with lock:
            with h5py.File(out, 'a') as h5:
                rep_group = h5[target_cond][sample_id]
                if rep_group.get(ref_name) is None:
                    rep_group.create_group(ref_name)
                ref_group = rep_group[ref_name]
                ref_group['forward'] = forward_comp
                ref_group['reverse'] = reverse_comp

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('scidx_fn', help='scidx file')
    parser.add_argument('bed_fns', help='file containing list of bed files')
    parser.add_argument('chrom_sizes', help='chrom.sizes file')
    parser.add_argument('out', help='output h5 file')
    parser.add_argument('--control', '-c', action='store_true', help='output control pileup')
    parser.add_argument('--ref_pattern', '-rp', default='(.+)_(\d+)bp.bed',
                        help='regular expression pattern to extract the reference name row width from the bed file name')
    parser.add_argument('--sample_pattern', '-sp', default='(\d+)_(.+?)_i5006_BY4741_-_YPD_(.+?)_XO_FilteredBAM',
                        help='regular expression pattern to extract the ID, target, and condition from the scidx file name')
    args = parser.parse_args()

    with open(args.bed_fns) as f:
        bed_fns = f.read().strip().split('\n')
    pileup(args.scidx_fn, bed_fns, args.chrom_sizes, args.out, control=args.control,
           ref_pattern=re.compile(args.ref_pattern), sample_pattern=re.compile(args.sample_pattern))