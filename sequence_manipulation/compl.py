compl_ = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
          'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
def compl(seq, reverse=True):
    return ''.join(map(compl_.__getitem__, reversed(seq) if reverse else seq))
