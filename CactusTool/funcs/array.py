"""
Dataset
"""
import numpy as np

def arithmetic_progression(a):
    """
    Arithmetic Progressions

    :param list a: Arithmetic Progressions
    :return: a0 and da
    """
    a = sorted(a)
    da = set(np.diff(a))
    if len(da) != 1:
        raise ImportError("{} is not a arithmetic progressions".format(a))
    return a[0], a[-1], da

def is_monotonically_increasing(t):
    """
    Wether the data is strictly increasing.
    
    :param t: Data samples
    :rtype:   bool
    """
    assert len(t) > 1, "A list has at least two Numbers"
    return np.diff(t).min() > 0

def list_bounds(x):
    x = np.asarray(x) 
    dx = x[1] - x[0]
    n = len(x)
    return np.linspace(x[0]-dx/2., x[-1]+dx/2, n+1)