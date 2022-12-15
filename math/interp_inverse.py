import numpy as np
import scipy.interpolate


def invert(f, x, kind='linear', vectorized=False):
    """
    Invert a function numerically
    Args:
        f: Function to invert
        x: Domain to invert the function on
        kind: Specifies the kind of interpolation as a string ('linear', 'spline', 'nearest', 'zero', 'slinear',
            'quadratic', 'cubic', 'previous', 'next')
        vectorized: Specifies if the input function is vectorized

    Returns:
        Inverted function with attributes x_min and x_max representing the domain bounds on which the function is defined
    """
    x = np.array(x)
    if not np.issubdtype(x.dtype, np.number):
        raise ValueError('Input domain is not numeric')

    if vectorized:
        y = f(x)
    else:
        y = np.array([f(x_i) for x_i in x])

    if not np.issubdtype(y.dtype, np.number):
        raise ValueError('Input function is not numeric')

    diff_signs = np.sign(np.diff(y))
    if not np.all(diff_signs == 1) and not np.all(diff_signs == -1):
        raise ValueError('Non-invertible function')

    if kind == 'spline':
        tck = scipy.interpolate.splrep(y, x)

        def f_inverse(x_new):
            return scipy.interpolate.splev(x_new, tck)

    else:
        f_inverse = scipy.interpolate.interp1d(y, x, kind=kind)

    f_inverse.x_min = min(y)
    f_inverse.x_max = max(y)
    return f_inverse


def invertnd(f, x, *other_vars, kind='linear', vectorized=False):
    """
    Invert a multivariate function numerically
    Args:
        f: Function to invert
        x: Domain to invert the function on (range of inverted function)
        *other_vars: Domain to invert the function on (parameters of inverted function)
        kind: Specifies the kind of interpolation as a string ('linear', 'nearest', 'cubic')
            (cubic only available for 1 or 2 variables)
        vectorized: Specifies if the input function is vectorized

    Returns:
        Inverted function where the first argument corresponds to the output of the original function
    """
    n = len(x)

    reshape_dim = np.ones(len(other_vars) + 1, dtype=int)
    reshape_dim[0] = n
    x_reshape = np.reshape(x, reshape_dim)
    reshape_dim[0] = 1

    if not np.issubdtype(x_reshape.dtype, np.number):
        raise ValueError('Input domain is not numeric')

    dim = [1, *(len(v) for v in other_vars)]
    x_arr = np.tile(x_reshape, dim)
    dim[0] = n

    v_arrs = []
    for i, v in enumerate(other_vars):
        reshape_dim[i + 1] = len(v)
        v_reshape = np.reshape(v, reshape_dim)
        reshape_dim[i + 1] = 1

        if not np.issubdtype(v_reshape.dtype, np.number):
            raise ValueError('Input domain is not numeric')

        dim[i + 1] = 1
        v_arrs.append(np.tile(v_reshape, dim))
        dim[i + 1] = len(v)

    if vectorized:
        y = f(x_arr, *v_arrs)
    else:
        def recursive_f(x_in, *v_in):
            if hasattr(x_in, '__iter__'):
                return [recursive_f(x_n, *v_n) for x_n, v_n in zip(x_in, zip(*v_in))]
            return f(x_in, *v_in)
        y = np.array(recursive_f(x_arr, *v_arrs))

    if not np.issubdtype(y.dtype, np.number):
        raise ValueError('Input function is not numeric')

    points = np.array(list(zip(y.flat, *(v.flat for v in v_arrs))))
    values = np.array(x_arr.flat)

    def f_inverse(x_new, *v_new):
        return scipy.interpolate.griddata(points, values, (x_new, *v_new), method=kind)

    return f_inverse
