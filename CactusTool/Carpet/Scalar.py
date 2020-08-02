"""
Cactus dataset main produced by `Carpet <https://carpetcode.org>`_, which is an adaptive mesh refinement and multi-patch driver for the Cactus. This modular processes output of CarpetIOScalar.
"""

from ..funcs import columns_asc
# from ..Lib.pygwanalysis import DiscreteFunction
import pandas as pd
import numpy as np
import os
import re


class CarpetIOScalar:
    """
    Thorn CarpetIOScalar provides I/O methods for outputting scalar values in ASCII format into files. This class can handle all CarpetIOScalar output. Read the dataset from files when given specified reduction operation.

    :param list files: A list of file in absolute path.
    :param str kind: Type of reduction operations.
    """
    def __init__(self, files):
        self.fname = {}
        pat = re.compile('^([a-zA-Z0-9_-]+)\.(minimum|maximum|norm1|norm2|average)?\.asc$')
        for file in files:
            mp = pat.match(os.path.basename(file))
            if mp is not None:
                thorn_var = mp.group(1)
                self.fname.update({thorn_var: file})

    def __getitem__(self, key):
        assert key in self.fname.keys(), "{} is not exist".format(key)
        return Variable(self.fname[key])

    def __contains__(self, key):
        return key in self.fname.keys()

class Variable:
    """
    For scalar, we don't need consider grid structure. Sometimes you want to process a vector or a tensor. :py:class:`Variable` can also handel it.

    :param list varfiles: A dict about variable and its file, this variable may be a vector or tensor.
    """
    def __init__(self, file):
        column = columns_asc(file)
        self.vars = column[2:]
        data = np.loadtxt(file, comments="#")
        self.dset = pd.DataFrame(data, columns=column)
        
    def Preview(self, var=None):
        dset = self.dset.drop('iteration', axis=1)
        if var:
            dset.plot(x='time', y=var)
        else:
            dset.plot(x='time', subplots=True)