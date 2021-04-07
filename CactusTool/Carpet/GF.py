"""
Cactus dataset main produced by `Carpet <https://carpetcode.org>`_, which is an adaptive mesh refinement and multi-patch driver for the Cactus. This modular processes output of CarpetIOHDF5 or CarpetIOASCII.
"""
from ..utils.file import header_h5, read, columns_asc
from ..utils.json import add_two_key
# from ..utils.log import logger
import numpy as np
import re
import os

class CarpetGF:
    """
    A class representing grid function.
    """
    def __init__(self, files, dim, ftype):
        self.files = files
        self.dim = dim
        self.ftype = ftype
        self.fname = set()
        pat = re.compile('^([a-zA-Z0-9_-]+)(\.\d+)?\.([xyz]*)(\.file_\d+)?\.(asc|h5)$')
        for f in files:
            mp = pat.match(os.path.basename(f))
            if mp is not None:
                thorn_var = mp.group(1)
                self.fname.add(thorn_var)
    
    def __getitem__(self, key):
        """
        Read the dataset when given specified file name.

        :param str key: file name.
        """
        assert key in self.fname, "{} is not exist".format(key)
        fileList = [f for f in self.files if os.path.basename(f).split(".")[0] == key]
        if self.ftype == 'h5':
            return CarpetIOHDF5(fileList, self.dim)
        elif self.ftype == 'asc':
            assert len(self.dim) == 1, "We do not recommend using ASCII to draw high-dimensional graphs"
            return CarpetIOASCII(fileList, self.dim)

    def __contains__(self, key):
        return key in self.fname


class CarpetIOHDF5:
    """
    Thorn CarpetIOHDF5 provides I/O methods for outputting gird function in HDF5 format into files. This class can handle all CarpetIOHDF5 output. Read the dataset from files and return all useful information stored in these files. We recommend use CarpetIOHDF5 to output all 2-D and 3-D grid function.
    """
    
    def __init__(self, files, dim):
        """
        :param list files: can either be a list of filenames or a single 
        :param str dim: dim
        """
        self.dim = dim
        self.vars = set()
        self.it = set()
        self.hierarchy = {}
        self.header = {}
        for file in files:
            self.header[file] = header_h5(file)
            for item in self.header[file].keys():
                p = self.header[file][item]
                self.vars.add(p['varname'])
                self.it.add(p['iteration'])
                rl = p['rl']
                if rl not in self.hierarchy:
                    self.hierarchy[rl] = set()
                if 'c' in p:
                    self.hierarchy[rl].add(p['c'])


    def dsets(self, var, it=-1, rl=-1, c=-1):
        """
        Select specified datasets.

        :param int it: -1: all iteration number.
        :param int rl: -1: all refinement level.
        :param int c: -1: all component.
        :return: a dict of datasets
        """
        assert var in self.vars, "Var {} is not exist".format(var)
        if it != -1:
            assert it in self.it, "iteration {} is not exist".format(it)
        if rl != -1:
            assert rl in self.hierarchy, "refinement level {} is not exist".format(rl)
        if c != -1:
            assert rl != -1
            assert c in self.hierarchy[rl], "component {} is not exist in refinement level {}".format(c, rl)

        p = {}
        for file in self.header:
            with read(file) as f: 
                for header in self.header[file]:
                    item = self.header[file][header]
                    if item['varname'] != var:
                        continue
                    if it != -1 and item['iteration'] != it:
                        continue
                    if rl != -1 and item['rl'] != rl:
                        continue
                    if c != -1 and item['c'] != c:
                        continue

                    mesh = f[header]
                    time = round(mesh.attrs.get('time', None), 6)
                    if time not in p:
                        p[time] = {}
                    rlevel = item['rl']
                    if rlevel not in p[time]:
                        p[time][rlevel] = {}
                    if 'c' in item:
                        component = item['c']
                    else:
                        component = 0
                    if component not in p[time][rlevel]:
                         p[time][rlevel][component] = {}
                    p[time][rlevel][component].update({'origin': mesh.attrs.get('origin', None)})
                    p[time][rlevel][component].update({'delta': mesh.attrs.get('delta', None)})
                    p[time][rlevel][component].update({'ghostzones': mesh.attrs.get('ghostzones', None)})
                    p[time][rlevel][component].update({'data': np.array(mesh)})

        return p


# class CarpetIOASCII:
#     """
#     Thorn CarpetIOASCII provides I/O methods for outputting gird function in ASCII format into files. This class can handle all CarpetIOASCII output. Read the dataset from files and return all useful information stored in these files.

#     :param str files: a single filename
#     """
#     def __init__(self, file, dim):
#         self.dim = dim
#         self.columns = columns_asc(file)
#         self.vars = self.columns[12:]
#         self._dsets = np.loadtxt(file, comments="#")

#     @property
#     def it(self):
#         it = self._dsets.T[0].astype(np.int)
#         return sorted(list(set(it)))

#     def dsets(self, var='data', it=0):
#         p = {}
#         dsets = self._dsets[self._dsets[:, 0] == it]
#         for rl in set(dsets[:, 2]):
#             dset = dsets[dsets[:, 2] == rl]
#             for c in set(dset[:, 3]):
#                 data = dset[dset[:, 3] == c].T
#                 p_tem = {}
#                 if self.vars[-1] == 'data':
#                     p_tem['data'] = data[12]
#                 else:
#                     assert var in self.vars, "Pick one from {}, because one file per group.".format(self.vars)
#                     p_tem['data'] = data[self.columns.index(var)]
#                 p_tem['coord'] = data[self.columns.index(self.dim)]
#                 p_tem['time'] = data[8][0]
#                 add_two_key(p, rl, c, p_tem)

#         return p
