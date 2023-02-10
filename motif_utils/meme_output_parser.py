import pandas as pd
import re

def parse_meme_output(meme_fn):
    """
    Parse MEME output file
    Args:
        meme_fn: Path to MEME output file

    Returns:
        List of DataFrames, one for each motif
    """

    wpattern = re.compile('MOTIF \w+ MEME-\d+\s+width =\s+(\d+)')
    hpattern = re.compile('Motif \w+ MEME-\d+ sites sorted by position p-value')
    datapattern = re.compile('(chr[\dM]+):(\d+)-(\d+)\(\.\)\s+([\+\-])\s+(\d+)\s+([\d\.\-e]+)\s+[\.\w]+\s+(\w+)')
    meme_dicts = []
    with open(meme_fn, 'r') as f:
        n = 0
        for line in f:
            match = re.search(wpattern, line)
            if match:
                n += 1
                meme_dicts.append({'chromosome': [], 'start': [], 'stop': [], 'strand': [], 'p-value': [], 'seq': []})
                width = int(match.group(1))
            if re.search(hpattern, line):
                f.readline()
                f.readline()
                f.readline()
                line = f.readline()
                while not line.startswith('-'):
                    chrom, start, stop, strand, idx, pval, seq = re.search(datapattern, line).groups()
                    start = int(start)
                    stop = int(start)
                    idx = int(idx)
                    pval = float(pval)
                    stop = start + idx + width - 1
                    start = start + idx - 1
                    meme_dicts[-1]['chromosome'].append(chrom)
                    meme_dicts[-1]['start'].append(start)
                    meme_dicts[-1]['stop'].append(stop)
                    meme_dicts[-1]['strand'].append(strand)
                    meme_dicts[-1]['p-value'].append(pval)
                    meme_dicts[-1]['seq'].append(seq)
                    line = f.readline()

    return [pd.DataFrame(meme_dict) for meme_dict in meme_dicts]
