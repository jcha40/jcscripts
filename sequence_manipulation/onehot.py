import numpy as np

def onehot(seq, alphadict={'A': np.array([1, 0, 0, 0]), 'C': np.array([0, 1, 0, 0]), 'G': np.array([0, 0, 1, 0]),
                           'T': np.array([0, 0, 0, 1]), 'N': np.array([0, 0, 0, 0])}):
    return np.array([alphadict[letter] for letter in seq])
