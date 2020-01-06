"""
The :py:mod:`~.outputfile` module module provides functions to load data and header in Cactus formats.
"""

from .dictionary import *
import numpy as np
import h5py
import bz2
import gzip
import re


def read(file):
    """
    Use different way to open data file.

    :Usage:

        >>> with read(file) as f:
    """
    if file.endswith('.bz2'):
        f = bz2.BZ2File(file)
    elif file.endswith('.gz'):
        f = gzip.GzipFile(file)
    elif file.endswith('.h5'):
        f = h5py.File(file, 'r')
    else:
        f = open(file)
    return f

def HDF5_header(file):
    """
    Return a dictionary about header information
    """
    pattern = re.compile(r'([^:]+)::(\S+) it=(\d+) tl=(\d+)( m=(\d+))? rl=(\d+)( c=(\d+))?')
    p = {}
    with read(file) as f:
        for item in list(f):
            m = pattern.match(item)
            if m == None:
                continue
            p[item] = {'thorn': m.group(1), 
                       'varname': m.group(2),
                       'iteration': int(m.group(3)),
                       'timelevel': int(m.group(4)),
                       'rl': int(m.group(7))}
            if m.group(6) != None:
                p[item].update({'ml': int(m.group(6))})
            if m.group(9) != None:
                p[item].update({'c': int(m.group(9))})

    return p

def HDF5_gf(file, slice):
    """
    Fetch data in .h5 file Given slice
    """
    data = dict()
    if 'NaNmask' in file:
        plane = 'xyz'
    else:
        try:
            plane, = re.match("\S*\.([xyz]*)\.h5", file).groups()
        except:
            RuntimeError("filename don't have plane information: %s" % (file))
    
    with read(file) as f:
        if slice not in sorted(list(f)):
            raise RuntimeError("%s not in %s" % (slice, file))

        dataset = f[slice]
        # for item in list(dataset.attrs):
        #     data[item] = dataset.attrs[item]
        data['plane'] = plane
        delta  = dataset.attrs.get('delta', None)
        origin = dataset.attrs.get('origin', None)
        size   = dataset.shape
        dim = len(size)
        for i in range(dim):
            data[plane[i]] = np.arange(0, size[(dim-1)-i])*delta[i] + origin[i]
        data['varname'] = str(dataset.attrs.get('name', None)).split("'")[1]
        data['it'] = dataset.attrs.get('timestep', None)
        data['rl'] = dataset.attrs.get('level', None)
        data['time'] = dataset.attrs.get('time', None)
        data['nghostzones'] = dataset.attrs.get('cctk_nghostzones', None)
        data['data'] = np.array(dataset)  
    
    return data

