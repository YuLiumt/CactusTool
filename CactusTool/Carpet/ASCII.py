"""
Cactus dataset main produced by `Carpet <https://carpetcode.org>`_, which is an adaptive mesh refinement and multi-patch driver for the Cactus. This modular processes output of CarpetIOASCII.
"""

from ..funcs import *
import pandas as pd
import numpy as np
import os
import re


class CarpetIOASCII:
    """
    Thorn CarpetIOASCII provides I/O methods for outputting gird function in ASCII format into files. This class can handle all CarpetIOASCII output. Read the dataset from files and return all useful information stored in these files.

    :py:class:`CarpetIOASCII` itself do nothing. You need specify the dimensions of gird function. This will pointer to :py:class:`Griddim`.

    * :py:attr:`CarpetIOASCII.x` x dimension of the gird function
    * :py:attr:`CarpetIOASCII.y` y dimension of the gird function
    * :py:attr:`CarpetIOASCII.z` z dimension of the gird function
    * :py:attr:`CarpetIOASCII.xy` xy dimension of the gird function
    * :py:attr:`CarpetIOASCII.xz` xz dimension of the gird function
    * :py:attr:`CarpetIOASCII.yz` yz dimension of the gird function
    * :py:attr:`CarpetIOASCII.xyz` xyz dimension of the gird function

    :param list files: A list of file in absolute path.
    """
    def __init__(self, files):
        self.files = ensure_list(files)
        self.x     = Griddim(self.files, 'x')
        self.y     = Griddim(self.files, 'y')
        self.z     = Griddim(self.files, 'z')
        self.xy    = Griddim(self.files, 'xy')
        self.xz    = Griddim(self.files, 'xz')
        self.yz    = Griddim(self.files, 'yz')
        self.xyz   = Griddim(self.files, 'xyz')

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
        self.files = [file for file in files if pat_fn.match(file).group(1) == self.dim]

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

    def __getitem__(self, key):
        p = {}
        for k, v in self.vars.items():
            if key in k:
                p[k] = v
        if len(p) > 1:
            print("Please make sure %s belong the same group" % (p.keys()))
        elif len(p) == 0:
            raise Exception("{} is not exist in reduction {}".format(key, self.kind))
        return Variable(p, self.dim)
    
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
    def __init__(self, varfiles, dim):
        self.varfiles = varfiles
        self.dim = dim
        self.Table = None

    def dataset(self, key=True, **kwargs):
        """
        CarpetIOASCII Dataset will store in pandas dataframe. Because the ASCII data structure more like a table.
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

    def grid_hierarchies(self, time=0.0):
        """
        Describes the geometry of the refined grid hierarchies, such as component number, ghost zones and refinement level. These all get from :py:meth:`Variable.dataset`. Grid hierarchies may change in the evolution.
        
        # So you need specify the time, the default is initial grid hierarchies (:math:`t = 0`).

        :param float time: grid hierarchies at time
        """
        assert time in var_ascii.Table['time'].values, "No such time value: {}".format(time)
        table = self.Table[self.Table['time'] == time]

        p = {}
        p['rl'] = self.Table['rl'].unique().astype(int)
        for i in self.dim:
            coord = self.Table[i]
        # self.c = self.Table['c'].unique().astype(int)
        #         column_header = columns(file)
        # table = pd.DataFrame(data, columns=column_header)
        # # check following column is what you want
        # column = ['it', 'rl', 'c', 'time'] + list(self.dim) + list(self.varfiles.keys())
        # p.append(table[column])

        return p

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
        output = "Time: [{}, {}]\n".format(min['time'], max['time'])
        output += "refinement level: {}\n".format(self.rl)
        output += "component: {}\n".format(self.c)
        for column in columns[12:]:
            output += "{}: Min is {:.4e};\tMean is {:.4e};\tMax is {:.4e}\n".format(column, min[column], mean[column], max[column])
        return output

def columns(file):
    """
    Fetch CarpetIOASCII header information.

    :param str file: file in absolute path
    :return: The column in given file.
    """
    with read(file) as f:
        columns=[]
        for line in f.readlines():
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
