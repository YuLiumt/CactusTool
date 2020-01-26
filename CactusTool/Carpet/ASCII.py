"""
Cactus dataset main produced by `Carpet <https://carpetcode.org>`_, which is an adaptive mesh refinement and multi-patch driver for the Cactus. This modular processes output of CarpetIOASCII.
"""

from ..funcs import *
from .grid import AMRGrid
import pandas as pd
import numpy as np
import os
import re


class CarpetIOASCII:
    """
    Thorn CarpetIOASCII provides I/O methods for outputting gird function in ASCII format into files. This class can handle all CarpetIOASCII output. Read the dataset from files and return all useful information stored in these files.

    :py:class:`CarpetIOASCII` itself do nothing. You need specify the dimensions of gird function. This will pointer to :py:class:`Griddim`.

    :param list files: can either be a list of filenames or a single filename
    """
    def __init__(self, files):
        self.files = ensure_list(files)
        self.dims = ['x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz']

    @property
    def available_dims(self):
        """
        available dimension in this dataset.

        :return: list of available dimension
        """
        p = []
        for dim in self.dims:
            if str(self[dim]):
                p.append(dim)
        return p

    def __getattr__(self, dim):
        """
        Specify the dimension of the gird function

        :param str dim: the dimension
        """
        assert dim in self.dims, "Does not include {} dim on ASCII".format(dim)
        return Griddim(self.files, dim)

    def __getitem__(self, dim):
        """
        Specify the dimension of the gird function by item

        :param str dim: the dimension
        """
        assert dim in self.dims, "Does not include {} dim on ASCII".format(dim)
        return Griddim(self.files, dim)

    def __str__(self):
        return "%s%s%s%s%s%s%s" % (self.x, self.y, self.z, self.xy, self.xz, self.yz, self.xyz)

class Griddim:
    """
    For gird function, We need choose the dimension of the gird function. Then choose the variable we want to processes, this will pointer to :py:class:`Variable`.

    :param list files: A list of file in absolute path.
    :param str dim: dimension
    """
    def __init__(self, files, dim):
        self.dim = dim
        pat_fn = re.compile("\S*\.([xyz]*)\.asc(\.(gz|bz2))?$")
        self.files = []
        for file in files:
            try: 
                if pat_fn.match(file).group(1) == self.dim:
                    self.files.append(file)
            except AttributeError:
                continue

    @property
    def vars(self):
        """
        All available variable in such dim.

        :return: A dict. key is the variable name, value is corresponding files.
        """
        Vars = {}
        for file in self.files:        
            with read(file) as f:
                vars=[]
                for line in f.readlines():
                    if "# data columns: " in line:
                        vars = vars + line.split()[3:]
                        break
            assert len(vars) > 0, "{}'s header fail to identify.".format(os.path.basename(file))

            for c, name in enumerate(vars):
                Vars.setdefault(name.split(":")[1], []).append(file)

        return Vars

    @property
    def available_variables(self):
        """
        All available variables in a given dimension.

        :return: list of available variables
        """
        return list(self.vars.keys())

    def __getitem__(self, key):
        assert key in self.available_variables, "{} is not exist in reduction {}".format(key, self.dim)
        return Variable(key, self.vars[key], self.dim)
    
    def __contains__(self, key):
        return key in self.vars

    def __str__(self):
        if self.vars:
            return "Available grid function with %s dimension:\n%s\n" % (str(self.dim).lower(), list(self.vars.keys()))
        else:
            return ""


class Variable:
    """
    In most case, the variable data have refined grid hirachies and may store in different files. We need combine them as needed. Sometimes you want to process a vector or a tensor. :py:class:`Variable` can also handel it.

    :param list varfiles: A dict about variable and its file, this variable may be a vector or tensor.
    :param str dim: dimension

    * :py:attr:`Variable.Table` source data
    """   
    def __init__(self, var, files, dim):
        self.var = var
        self.files = files
        self.dim = dim

    @property
    def dataset(self):
        """
        CarpetIOASCII Dataset will store in pandas dataframe. Because the ASCII data structure more like a table. This dataset will store in :py:attr:`Variable.dataset`.

        :return: DataFrame
        """
        tem = []
        for file in self.files:
            data = np.loadtxt(file, comments="#")
            column = columns_asc(file)
            tem_c = pd.DataFrame(data, columns=column)
            tem.append(tem_c)
        # concat different component in same variable
        return pd.concat(tem).drop_duplicates()

    @property
    def it(self):
        return self.dataset['it'].unique().astype(int).tolist()

    @property
    def time(self):
        return pd.Series(self.dataset['time'].values, index=self.dataset['it'].astype(int)).drop_duplicates()
        
    def temporary(self, it=0):
        dataset = self.dataset[self.dataset.it == it]
        column = ['rl', 'c', self.var]
        for dim in self.dim:
            column += dim
        dset = pd.DataFrame(dataset, columns=column) 
        return AMRGrid(dset, self.dim, self.var, 'ascii')

    def value(self, time=0):
        dataset = self.dataset[self.dataset.time == time]
        points = tuple([dataset[dim].values for dim in self.dim])
        value = dataset[self.var].values
        return LinearNDInterpolator(points, value)
        # return points


    