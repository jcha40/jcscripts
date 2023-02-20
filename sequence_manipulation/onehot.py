import numpy as np

def onehot(seq, alphadict={'A': (1, 0, 0, 0), 'C': (0, 1, 0, 0), 'G': (0, 0, 1, 0), 'T': (0, 0, 0, 1), 'N': (0, 0, 0, 0)}):
    return np.array([alphadict[letter] for letter in seq])
