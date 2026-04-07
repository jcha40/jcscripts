from numpy import exp, log, log1p, ndarray, empty_like, float64
from scipy.special import betaln

_loghalf = log(.5)
def logbetah(a, b, n=20):
    '''
    Compute the logarithm of the regularized half-beta function using the continued fraction representation.
    '''
    if not isinstance(a, (int, float, ndarray)) or not isinstance(b, (int, float, ndarray)):
        raise TypeError('a and b must be int, float, or ndarray')

    if isinstance(a, ndarray) and isinstance(b, ndarray) and a.shape != b.shape:
        raise ValueError('a and b must have the same shape')

    if not isinstance(n, int) or n <= 0:
        raise ValueError('n must be a positive integer')
    
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if a < b:
            return log1p(-exp(logbetah(b, a, n=n)))
        frac = 1.
        for k in range(n, 0, -1):
            frac = 1 - (((a + k - 1) * (a + b + k - 1) * .5) / ((a + k * 2 - 2) * (a + k * 2 - 1))) / (1 + ((k * (b - k) * .5) / ((a + k * 2 - 1) * (a + k * 2))) / frac)
        return a * _loghalf + b * _loghalf - log(a) - betaln(a, b) - log(frac)

    ans = empty_like(a + b, dtype=float64)
    mask = a < b
    if isinstance(a, (int, float)):
        a_lt = a_ge = a
    else:
        a_lt = a[mask]
        a_ge = a[~mask]
    if isinstance(b, (int, float)):
        b_lt = b_ge = b
    else:
        b_lt = b[mask]
        b_ge = b[~mask]
    
    if mask.any():
        ans[mask] = log1p(-exp(logbetah(b_lt, a_lt, n=n)))

    frac = 1.
    for k in range(n, 0, -1):
        frac = 1 - (((a_ge + k - 1) * (a_ge + b_ge + k - 1) * .5) / ((a_ge + k * 2 - 2) * (a_ge + k * 2 - 1))) / (1 + ((k * (b_ge - k) * .5) / ((a_ge + k * 2 - 1) * (a_ge + k * 2))) / frac)
    ans[~mask] = a_ge * _loghalf + b_ge * _loghalf - log(a_ge) - betaln(a_ge, b_ge) - log(frac)

    return ans