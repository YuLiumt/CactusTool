"""
Cactus dataset main produced by `Carpet <https://carpetcode.org>`_, which is an adaptive mesh refinement and multi-patch driver for the Cactus. This modular processes output of CarpetIOScalar.
"""

from ..funcs import *
import matplotlib.pyplot as plt
import numpy as np
import os
import re


class CarpetIOScalar:
    """
    Thorn CarpetIOScalar provides I/O methods for outputting scalar values in ASCII format into files. This class can handle all CarpetIOScalar output. Read the dataset from files and return all executable reduction operation.

    :py:class:`CarpetIOScalar` itself do nothing. You need specify the reduction operation. This will pointer to :py:class:`ScalarReduction`.

    * :py:attr:`CarpetIOScalar.min` The minimum of the gird function
    * :py:attr:`CarpetIOScalar.max` The maximum of the gird function
    * :py:attr:`CarpetIOScalar.norm1` The norm1 of the gird function
    * :py:attr:`CarpetIOScalar.norm2` The norm2 of the gird function
    * :py:attr:`CarpetIOScalar.average` The average of the gird function
    * :py:attr:`CarpetIOScalar.none` The variable is scalar, don't need do any reduction

    :param list files: A list of file in absolute path.

    >>> ds = CarpetIOScalar(Scalar_file)
    >>> print(ds)
    Available min timeseries:
    []
    Available max timeseries:
    ['H']
    """
    def __init__(self, files):
        self.files    = ensure_list(files)
        self.min      = ScalarReduction(self.files, 'minimum')
        self.max      = ScalarReduction(self.files, 'maximum')
        self.norm1    = ScalarReduction(self.files, 'norm1')
        self.norm2    = ScalarReduction(self.files, 'norm2')
        self.average  = ScalarReduction(self.files, 'average')
        self.none     = ScalarReduction(self.files, None)

    def __str__(self):
        return "%s%s%s%s%s%s%s" % (self.path, self.scalar, self.min, self.max, self.norm1, self.norm2, self.average)


class ScalarReduction:
    """
    For gird function, We need choose which type of reduction operations.

    :param list files: A list of file in absolute path.
    :param str kind: Type of reduction operations.
    """
    def __init__(self, files, kind):
        pat_fn = re.compile("\S*\.(minimum|maximum|norm1|norm2|average)?\.asc(\.(gz|bz2))?$")
        self.kind = kind
        self.files = [file for file in files if pat_fn.match(file).group(1) == self.kind]

    @property
    def vars(self):
        """
        All available variable in such reduction.

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

    def __getitem__(self, key):
        vars = [var for var in self.vars if key in var]
        if len(vars) > 1:
            print("Please make sure %s belong the same group" % (vars))
        elif len(vars) == 0:
            raise Exception("{} is not exist in reduction {}".format(key, self.kind))
        files = []
        for var in vars:
            files += self.vars[var]
        return Variable(files, vars)
    
    def __contains__(self, key):
        return key in self.vars
        
    def __str__(self):
        if self.vars:
            return "Available %s timeseries:\n%s\n" % (str(self.kind).lower(), list(self.vars.keys()))
        else:
            return None


class Variable:
    """
    In most case, the variable data have refined grid hirachies and may store in different files. We need combine them as needed. Sometimes you want to process a vector or a tensor. :py:class:`Variable` can also handel it.

    :param list files: A list of file in absolute path.
    :param list var: The variable you want to deal with, may be a vector or tensor.
    """
    def __init__(self, files, vars):
        self.files = ensure_list(files)
        self.var = ensure_list(vars)


        self.iteration = {}
        for file in self.files:
            self.iteration.setdefault(iteration(file), []).append(file)
        self.column = column_header(self.files[0])
        Alldata = merge_filedata(self.files)
        self.data = remove_duplicate_iters(Alldata)

    @property
    def t(self):
        i = self.column.index('time')
        return self.data[i, :]
    
    @property
    def y(self):
        i = self.column.index(self.var)
        return self.data[i, :]

    def __str__(self):
        return "{}".format(json.dumps(self.iteration, sort_keys=True, indent=4))


# def Data(file):
#     """
#     Get the data in Scalar file.

#     Args:
#         file (str): open file
    
#     Return:
#         data
#     """
#     return np.loadtxt(file, comments="#", unpack=True)

def remove_duplicate_iters(data):
    """
    Remove overlapping segments from a time series. According to the first column
    """
    u, indices = np.unique(data[0,:], return_index=True)
    return data[:, indices]

def merge_filedata(filelist):
    """
    Some Variable data locate in different components.

    Args:
        file (list): a list of scalar file with different components.
    
    Return:
        data contain all components.
    """
    # Get the data from the first file
    data = np.loadtxt(filelist[0], comments="#", unpack=True)
    for i in range(1, len(filelist)):
        try: 
            tmp = np.loadtxt(filelist[i], comments="#", unpack=True)
            data = np.append(data, tmp, axis=1)
        except:
            print("[ERROR] Unable to load file:", filelist[i])
    return data

def column_header(file):
    with read(file) as f:
        columns=[]
        for line in f.readlines():
            if "# 1:iteration 2:time 3:data" in line:
                columns = columns + line.split()[1:]
            if "# column format:" in line:
                columns = line.split()[3:]
            if "# data columns: " in line:
                del columns[-1]
                columns = columns + line.split()[3:]
                break
    if len(columns) > 0:
        return [name.split(":")[1] for c, name in enumerate(columns)]
    else:
        raise Exception("File: {} Header fail to identify.".format(file))

def iteration(file):
    scalar_pat = re.compile("\S*/(output-\d\d\d\d)/\S*\.(minimum|maximum|norm1|norm2|norm_inf|average)?\.asc(\.(gz|bz2))?$")
    iteration = scalar_pat.match(file)
    if iteration is not None:
        return iteration.group(1)
    else:
        return "output-0000"


