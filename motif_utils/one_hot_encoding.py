import numpy as np

def ohe_motif(motif, alphabet='ACGT'):
    """One-hot-encode a motif.

    Args:
        motif (str): Motif to encode.
        alphabet (str): Alphabet to use for encoding.

    Returns:
        np.array: One-hot-encoded motif.
    """
    motif = motif.upper()
    alphadict = {letter: i for i, letter in enumerate(alphabet)}
    ohe = np.zeros((len(motif), len(alphadict)), dtype=int)
    for i, letter in enumerate(motif):
        ohe[i, alphadict[letter]] = 1
    return ohe
