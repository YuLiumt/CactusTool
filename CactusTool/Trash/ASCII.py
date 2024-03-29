"""
Cactus dataset main produced by `Carpet <https://carpetcode.org>`_, which is an adaptive mesh refinement and multi-patch driver for the Cactus. This modular processes output of CarpetIOASCII.
"""

from ..utils import read, ensure_list, columns_asc
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
        pat_fn = re.compile("\S*\.([xyz]+)\.asc(\.(gz|bz2))?$")
        self.dims = set(pat_fn.match(file).group(1) for file in self.files)

    def __getattr__(self, dim):
        """
        Specify the dimension of the gird function

        :param str dim: the dimension
        """
        assert dim in self.dims, "Does not include {} dim on ASCII".format(dim)
        return self[dim]

    def __getitem__(self, dim):
        """
        Specify the dimension of the gird function by item

        :param str dim: the dimension
        """
        assert dim in self.dims, "Does not include {} dim on ASCII".format(dim)
        files = [file for file in self.files if ".{}.".format(dim) in file]
        return Griddim(files, dim)


class Griddim:
    """
    For gird function, We need choose the dimension of the gird function. Then choose the variable we want to processes, this will pointer to :py:class:`Variable`.

    :param list files: A list of file in absolute path.
    :param str dim: dimension
    """
    def __init__(self, files, dim):
        self.dim = dim
        self.files = files
        self.vars = {}

        result = map(self.var_info, self.files)
        for key, values in dict(zip(self.files, result)).items():
            for value in values:
                self.vars.setdefault(value, []).append(key)
        
    @staticmethod
    def var_info(file):
        vars = set()
        with read(file) as f:
            for line in f.readlines():
                if "# data columns: " in line:
                    for var in line.split()[3:]:
                        vars.add(var.split(":")[1]) 
                    break
            assert len(vars) > 0, "{}'s header fail to identify.".format(os.path.basename(file))

        return vars

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
        self.dataset = self._init_dataset()
 
    def _init_dataset(self):
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
        return sorted(self.dataset['it'].unique().astype(int).tolist())

    @property
    def time(self):
        return pd.Series(self.dataset['time'].values, index=self.dataset['it'].astype(int)).drop_duplicates().sort_index()
        
    def temporary(self, it=0):
        dataset = self.dataset[self.dataset.it == it]
        column = ['rl', 'c', self.var]
        for dim in self.dim:
            column += dim
        dset = pd.DataFrame(dataset, columns=column) 
        return AMRGrid(dset, self.dim, self.var)
