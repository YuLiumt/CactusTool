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
    da = set(np.diff(a))
    if len(da) != 1:
        raise ImportError("{} is not a arithmetic progressions".format(a))
    return da

def merge_data_simple(alldat):
    """Merges a list of RegData instances into one, assuming they
    all have the same grid spacing and filling a regular grid completely,
    with the same values in overlapping regions (ghost zones etc).
    Beware, the assumptions are not checked.
    """
    if len(alldat)==0:
        return None
    if len(alldat)==1:
        return alldat[0]

    mg    = merge_geom(alldat)
    data  = zeros(mg.shape(), dtype=alldat[0].data.dtype)

    for d in alldat:
        i0      = ((d.x0()-mg.x0())/mg.dx() + 0.5).astype(int32)
        i1      = i0 + d.shape()
        i       = [slice(j0,j1) for j0,j1 in zip(i0,i1)]
        data[i] = d.data
    #
    return RegData(mg.x0(),mg.dx(),data, reflevel=alldat[0].reflevel(), component=-1)

def CompData(a):
    """Composite data consisting of one or more regular datasets with 
    different grid spacings, i.e. a mesh refinement hirachy. The grid 
    spacings should differ by powers of two. Origins of the components 
    are shifted relative to each other only by multiples of the finest 
    spacing. Basic arithmetic operations are defined for this class, as 
    well as interpolation and resampling. This class can be iterated over 
    to get all the regular datasets, ordered by refinement level and
    componen number.
    """
    return a

