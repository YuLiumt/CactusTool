"""
The `cactus_scalars` module provides functions to load 
timeseries in Cactus formats and a class `ScalarsDir` for easy 
access to all timeseries in a Cactus simulation directory. This module 
is normally not used directly, but from the `simdir` module. 
The data loaded by this module is represented as `TimeSeries` objects.
"""

from .outputfile import read
import matplotlib.pyplot as plt
import numpy as np
import json
import os
import re


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


class Variable:
    def __init__(self, files, var):
        self.files = files
        self.var = var
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


def var_header(file):
    """
    All Scalar variable.
    """
    with read(file) as f:
        vars=[]
        for line in f.readlines():
            if "# data columns: " in line:
                vars = vars + line.split()[3:]
                break
    if len(vars) > 0:
        return [name.split(":")[1] for c, name in enumerate(vars)]
    else:
        raise Exception("File: {} Header fail to identify.".format(file))


class ScalarReduction:
    def __init__(self, files, kind):
        pat_fn = re.compile("\S*\.(minimum|maximum|norm1|norm2|average)?\.asc(\.(gz|bz2))?$")
        self.kind = kind
        self.files = [file for file in files if pat_fn.match(file).group(1) == self.kind]
        self.vars = {}
        for file in self.files:
            for var in var_header(file):
                self.vars.setdefault(var, []).append(file)

    def __getitem__(self, key):
        if key in self.vars:
            return Variable(self.vars[key], key)
        else:
            raise Exception("{} is not exist in reduction {}".format(key, self.kind))
    
    def __contains__(self, key):
        return key in self.vars
        
    def __str__(self):
        return "Available %s timeseries:\n%s\n" % (str(self.kind).lower(), list(self.vars.keys()))


class ScalarBase:
    """
    Operation in Carpet Scalar file.
    """
    def __init__(self, Sim):
        self.files    = Sim.scafiles
        self.min      = ScalarReduction(self.files, 'minimum')
        self.max      = ScalarReduction(self.files, 'maximum')
        self.norm1    = ScalarReduction(self.files, 'norm1')
        self.norm2    = ScalarReduction(self.files, 'norm2')
        self.average  = ScalarReduction(self.files, 'average')
        self.none     = ScalarReduction(self.files, None)

    def __str__(self):
        return "%s\n%s\n%s\n%s\n%s\n%s\n" % (self.min, self.max, self.norm1, self.norm2, self.average, self.none)