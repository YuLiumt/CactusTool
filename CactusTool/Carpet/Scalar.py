"""
Cactus dataset main produced by `Carpet <https://carpetcode.org>`_, which is an adaptive mesh refinement and multi-patch driver for the Cactus. This modular processes output of CarpetIOScalar.
"""

from ..funcs import *
import pandas as pd
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
        return "%s%s%s%s%s%s" % (self.min, self.max, self.norm1, self.norm2, self.average, self.none)


class ScalarReduction:
    """
    For gird function, We need choose which type of reduction operations. Then choose the variable we want to processes, this will pointer to :py:class:`Variable`.

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
        p = {}
        for k, v in self.vars.items():
            if key in k:
                p[k] = v
        if len(p) > 1:
            print("Please make sure %s belong the same group" % (p.keys()))
        elif len(p) == 0:
            raise Exception("{} is not exist in reduction {}".format(key, self.kind))
        return Variable(p)
    
    def __contains__(self, key):
        return key in self.vars
        
    def __str__(self):
        if self.vars:
            return "Available %s timeseries:\n%s\n" % (str(self.kind).lower(), list(self.vars.keys()))
        else:
            return ""


class Variable:
    """
    For scalar, we don't need consider grid structure. These data may store in different files, we just simple combine them and remove duplicate data. This will done by pandas DataFrame. Sometimes you want to process a vector or a tensor. :py:class:`Variable` can also handel it.

    :param list varfiles: A dict about variable and its file, this variable may be a vector or tensor.

    * :py:attr:`Variable.Table` source data
    """
    def __init__(self, varfiles):
        self.varfiles = varfiles
        self.Table = None

    def dataset(self, key=True, **kwargs):
        """
        CarpetIOScalar Dataset will store in pandas dataframe, because the ASCII data structure more like a table. This dataset will store in :py:attr:`Variable.Table`.

        :param kwargs: Unknown keyword arguments are passed to :py:func:`pd.concat()`.

        :return: DataFrame
        """
        files = []
        p = []
        for var in self.varfiles.keys():
            for file in self.varfiles[var]:
                filename = os.path.basename(file)
                if filename in files:
                    continue
                files.append(filename)
                data = np.loadtxt(file, comments="#")
                column = columns(file)
                p.append(pd.DataFrame(data, columns=column)) 

        if key:
            self.Table = pd.concat(p, keys=files)
        else:
            self.Table = pd.concat(p, **kwargs)

        return self.Table

    def preview(self, **kwargs):
        """
        :py:meth:`Variable.preview` just simple preview. We will use :py:func:`pandas.DataFrame.plot()` to do it. It is best to run :py:meth:`Variable.dataset` and check the dataset by eye before executing it.

        :param kwargs: Unknown keyword arguments are passed to :py:func:`pandas.DataFrame.plot()`.
        """
        if self.Table.empty:
            print("Use default method combine multi data file.")
            self.dataset()
        assert 'time' in self.Table, "Dataset don't have time column"
        self.Table.plot(x='time', **kwargs)

    def __str__(self):
        if self.Table.empty:
            self.dataset()

        columns = self.Table.columns
        min = self.Table.min()
        max = self.Table.max()
        mean = self.Table.mean()
        output = "Time: [{}, {}] with dt={}\n".format(min['time'], max['time'], self.Table['time'][1]-min['time'])
        for column in columns[2:]:
            output += "{}: Min is {:.4e};\tMean is {:.4e};\tMax is {:.4e}\n".format(column, min[column], mean[column], max[column])
        return output

def columns(file):
    """
    Fetch CarpetIOScalar header information.

    :param str file: file in absolute path
    :return: The column in given file.
    """
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
