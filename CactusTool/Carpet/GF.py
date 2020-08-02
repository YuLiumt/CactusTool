"""
Cactus dataset main produced by `Carpet <https://carpetcode.org>`_, which is an adaptive mesh refinement and multi-patch driver for the Cactus. This modular processes output of CarpetIOHDF5 or CarpetIOASCII.
"""
from ..funcs import header_h5, select_header_h5, read, add_two_key
import numpy as np
import re
import os

class CarpetGF:
    def __init__(self, files, dim, format):
        self.dim = dim
        self.format = format
        self.fname = {}
        pat = re.compile('^([a-zA-Z0-9_-]+)\.([xyz]*)(\.file_\d+)?\.(asc|h5)$')
        for file in files:
            mp = pat.match(os.path.basename(file))
            if mp is not None:
                thorn_var = mp.group(1)
                self.fname.setdefault(thorn_var, []).append(file)
    
    def __getitem__(self, key):
        assert key in self.fname.keys(), "{} is not exist".format(key)
        files = self.fname[key]
        if self.format == 'h5':
            return CarpetIOHDF5(files, self.dim)
        elif self.format == 'asc':
            raise Exception("CactusTool currently only support '.h5'!")

    def __contains__(self, key):
        return key in self.fname.keys()

class CarpetIOHDF5:
    """
    Thorn CarpetIOHDF5 provides I/O methods for outputting gird function in HDF5 format into files. This class can handle all CarpetIOHDF5 output. Read the dataset from files and return all useful information stored in these files. We recommend use CarpetIOHDF5 to output all 2-D and 3-D grid function.

    :param str file: a single filename
    """
    def __init__(self, files, dim):
        self.dim = dim
        self.header = {}
        for file in files:
            self.header[file] = header_h5(file)

    @property
    def vars(self):
        vars = set()
        for file in self.header:
            for item in self.header[file]:
                vars.add(self.header[file][item]['varname'])
        return list(vars)

    @property
    def it(self):
        iteration = set()
        for file in self.header:
            for item in self.header[file]:
                iteration.add(int(self.header[file][item]['iteration']))
        return sorted(list(iteration))

    def temporary(self, var, it=0):
        p = {}
        for file in self.header:
            with read(file) as f: 
                for header in select_header_h5(self.header[file], var, it=it):
                    rl = self.header[file][header]['rl']
                    c = self.header[file][header]['c']
                    dset = {}
                    mesh = f[header]
                    origin = mesh.attrs.get('origin', None)
                    delta = mesh.attrs.get('delta', None)
                    data = np.array(mesh)
                    size = mesh.shape
                    n = len(self.dim)
                    for i in range(n):
                        dset[self.dim[i]] = np.arange(0, size[(n-1)-i]+1)*delta[i] + origin[i] - delta[i]/2
                    dset['data'] = data
                    add_two_key(p, rl, c, dset)
        return p

#     @property
#     def rl(self):
#         rl = set()
#         for file in self.header:
#             for item in self.header[file]:
#                 rl.add(int(self.header[file][item]['rl']))
#         return sorted(list(rl))

#     def temporary(self, var, it=0, rl=-1):
#         dsets = []
#         if rl == -1:
#             rl = self.rl
#         for level in list(rl):
#             for file in self.header:
#                 with read(file) as f: 
#                     headers = select_header_h5(self.header[file], var, it=it, rl=level)
            
#                     for slice in headers:
#                         p = dict()
#                         mesh = f[slice]
#                         p['level'] = mesh.attrs.get('level', None)
#                         origin = mesh.attrs.get('origin', None)
#                         delta = mesh.attrs.get('delta', None)
#                         p['origin'] = origin
#                         data = np.array(mesh)
#                         size = mesh.shape
#                         n = len(self.dim)
#                         for i in range(n):
#                             p[self.dim[i]] = np.arange(0, size[(n-1)-i]+1)*delta[i] + origin[i] - delta[i]/2
#                         p['data'] = data
#                         dsets.append(p)
#         return AMRGrid(dsets, self.dim, var)

# class AMRGrid:
#     def __init__(self, dsets, dim, var):
#         self.dsets = dsets
#         self.dim = dim
#         self.var = var
    
#     def Preview(self, axlim=None, Normalize=None):
#         import matplotlib.pyplot as plt
#         import matplotlib.colors as colors

#         vmax = np.amax(self.dsets[0]['data'])
#         vmin = np.amin(self.dsets[0]['data'])
#         for item in self.dsets:
#             if vmax < np.amax(item['data']):
#                 vmax = np.amax(item['data'])
#             if vmin > np.amin(item['data']):
#                 vmin = np.amin(item['data'])
#         for item in self.dsets:
#             if Normalize is None:
#                 plt.pcolormesh(item[self.dim[0]], item[self.dim[1]], item['data'], vmin=vmin, vmax=vmax)
#             elif Normalize == 'LogNorm':
#                 assert vmin > 0, "values must all be positive"
#                 plt.pcolormesh(item[self.dim[0]], item[self.dim[1]], item['data'], norm=colors.LogNorm(vmin=vmin, vmax=vmax))

#         if axlim is not None:
#             plt.xlim(axlim)
#             plt.ylim(axlim)
#         plt.xlabel('x')
#         plt.ylabel('y')
#         plt.title(self.var)
#         plt.colorbar()
#         plt.show()