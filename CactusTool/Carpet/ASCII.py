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

    :param list files: can either be a list of filenames or a single filename
    """
    def __init__(self, files):
        self.files = ensure_list(files)

    @property
    def x(self):
        """
        x dimension of the gird function
        """
        return Griddim(self.files, 'x')

    @property
    def y(self):
        """
        y dimension of the gird function
        """
        return Griddim(self.files, 'y')

    @property
    def z(self):
        """
        z dimension of the gird function
        """
        return Griddim(self.files, 'z')

    @property
    def xy(self):
        """
        xy dimension of the gird function
        """
        return Griddim(self.files, 'xy')

    @property
    def xz(self):
        """
        xz dimension of the gird function
        """
        return Griddim(self.files, 'xz')

    @property
    def yz(self):
        """
        yz dimension of the gird function
        """
        return Griddim(self.files, 'yz')

    @property
    def xyz(self):
        """
        xyz dimension of the gird function
        """
        return Griddim(self.files, 'xyz')

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

    def __getitem__(self, key):
        p = {}
        for k, v in self.vars.items():
            if key in k:
                p[k] = v
        if len(p) == 0:
            raise Exception("{} is not exist in reduction {}".format(key, self.kind))
        return Variable(p)
    
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
    def __init__(self, varfiles):
        self.varfiles = varfiles
        self._init_dataset()

    def _init_dataset(self):
        """
        CarpetIOASCII Dataset will store in pandas dataframe. Because the ASCII data structure more like a table. This dataset will store in :py:attr:`Variable.dataset`.

        :return: DataFrame
        """
        files = []
        p = pd.DataFrame()
        for var in self.varfiles.keys():
            tem = []
            for file in self.varfiles[var]:
                filename = os.path.basename(file)
                if filename in files:
                    continue
                files.append(filename)
                data = np.loadtxt(file, comments="#")
                if data.size == 0:
                    continue
                column = columns_asc(file)
                tem_c = pd.DataFrame(data, columns=column)
                tem.append(tem_c)
            # concat different component in same variable
            if len(tem) == 0:
                continue
            else:
                tem_p = pd.concat(tem).drop_duplicates()
            # merge different variable
            if p.empty:
                p = tem_p
            else:
                p = pd.merge(p, tem_p, how='outer', on=['time','it'])
 
        self.dataset = p

    def grid_hierarchies(self, time=0.0):
        """
        # TODO:
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

    def slice(self, time, rl=-1, c=-1):
        """
        Employ slicing to select sets of data from :py:attr:`Variable.dataset`.

        :param float time: time
        :param int rl: refinement level. -1: all refinement level
        :param int c: component. -1: all component

        :return: slice of DataFrame
        """
        return None
