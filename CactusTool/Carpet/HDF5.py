"""
Cactus dataset main produced by `Carpet <https://carpetcode.org>`_, which is an adaptive mesh refinement and multi-patch driver for the Cactus. This modular processes output of CarpetIOHDF5.
"""

from ..funcs import header_h5, select_header_h5, dataset_yt, ensure_list, read
# from .grid import AMRGrid
# import multiprocessing as mp
import pandas as pd
import numpy as np
import os
import re
import yt

class CarpetIOHDF5:
    """
    Thorn CarpetIOHDF5 provides I/O methods for outputting gird function in HDF5 format into files. This class can handle all CarpetIOHDF5 output. Read the dataset from files and return all useful information stored in these files. We recommend use CarpetIOHDF5 to output all 2-D and 3-D grid function.

    :py:class:`CarpetIOHDF5` itself do nothing. You need specify the dimensions of gird function. This will pointer to :py:class:`Griddim`.

    :param list files: can either be a list of filenames or a single filename
    """
    def __init__(self, files):
        self.files = ensure_list(files)
        pat_fn = re.compile("\S*\.([xyz]*)(\.file_\d+)?\.h5(\.(gz|bz2))?$")
        self.dims = set(pat_fn.match(file).group(1) for file in self.files)

    def __getattr__(self, dim):
        """
        Specify the dimension of the gird function

        :param str dim: the dimension
        """
        assert dim in self.dims, "Grid function doesn't have dim {}".format(dim)
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
        parser = re.compile(r'([^:]+)::(\S+) it=(\d+) tl=(\d+)( m=0)? rl=(\d+)( c=(\d+))?')
        vars = set()
        with read(file) as f:
            for level in sorted(list(f)):
                m = parser.match(level)
                if m is not None:
                    vars.add(m.group(2)) 
        return vars

    @property
    def available_variables(self):
        """
        All available variables in a given dimension.

        :return: list of available variables
        """
        return list(self.vars.keys())

    def __getitem__(self, key):
        assert key in self.available_variables, "{} is not exist in dimention {}".format(key, self.dim)
        return Variable(key, self.vars[key], self.dim)

    def __contains__(self, key):
        return key in self.vars

    def __str__(self):
        if self.vars:
            return "Available grid function with %s dimension:\n%s\n" % (self.dim, list(self.vars.keys()))
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

        self.header = {}
        for file in self.files:
            header = header_h5(file)
            self.header.update(select_header_h5(header, self.var))

    @property
    def it(self):
        iteration = set()
        for item in self.header:
            iteration.add(int(self.header[item]['iteration']))
        return sorted(list(iteration))

    def temporary(self, it=0):
        headers = select_header_h5(self.header, self.var, it=it)
        dsets = []
        for item in headers:
            dsets.append(dataset_yt(headers[item]['file'], item, self.var))
        return yt.load_amr_grids(dsets, [1000, 1000, 1000])
        # return AMRGrid(dset, self.dim, self.var)
        # return dsets


    # @property
    # def it(self):
    #     return self.dataset['it'].unique().astype(int).tolist()

    # @property
    # def time(self):
    #     return pd.Series(self.dataset['time'].values, index=self.dataset['it'].astype(int)).drop_duplicates().sort_index()

    # @property
    # def time(self):
    #     p = []
    #     for item in self.dataset:
    #         time.append(self.dataset[item]['time'])
    #     return time

    # def temporary(self, it=0):
    #     headers = select_header_h5(self.dataset, self.var, it=it)
    #     dsets = []
    #     for item in headers:
    #         dsets.append(dataset_h5(self.dataset[item]['file'], item))
    #     # return AMRGrid(dset, self.dim, self.var, 'hdf5')
    #     return dsets


    # def grid_hierarchies(self): 
    #     """
    #     Describes the geometry of the refined grid hierarchies, such as component number, ghost zones and refinement level. Grid hierarchies may change in the evolution. These all get from the header of files.

    #     :return: a dict about grid_hierarchies
    #     """
    #     parser = re.compile(r'([^:]+)::(\S+) it=(\d+) tl=(\d+)( m=0)? rl=(\d+)( c=(\d+))?')
    #     for var in self.varfiles.keys():
    #         for file in self.varfiles[var]:
    #             filename = os.path.basename(file)
    #             if filename in files:
    #                 continue
    #     return None

    # def slice(self, meshgrid='HYDROBASE::press it=0 tl=0 rl=0 c=10'):
    #     """
    #     CarpetIOHDF5 is different with CarpetIOASCII. We don't need read all data in the beginning. 2-D or 3-D data is huge, reading all data at one time is a waste of resources.

    #     :return: DataFrame
    #     """
    #     with read(self.files[0]) as f:
    #         mesh = f[meshgrid]
    #         delta = mesh.attrs['delta']
    #         origin = mesh.attrs['origin']
    #         sizeA = mesh.shape
    #         tmpX = np.arange(0,sizeA[1])*delta[0]+origin[0]
    #         tmpY = np.arange(0,sizeA[0])*delta[1]+origin[1]

    #         grid = np.meshgrid(tmpX, tmpY)
    #         data = np.array(mesh) 
    #     return grid, data



# def merge_filedata(filelist):
#     p = []
#     for file in filelist:
#         with read(file) as f:
#             for dset in sorted(list(f)):
#                 infos = dict()
#                 REG = re.match('(\S+)::(\S+) it=(\d+)',dset)
#                 if REG:
#                     infos['group'] = REG.groups()[0]
#                     infos['var']   = REG.groups()[1]
#                     infos['it']    = int(REG.groups()[2])
#                 REG = re.search('tl=(\d+)',dset); 
#                 if REG: 
#                     infos['tl']=int(REG.groups()[0])
#                 REG = re.search('rl=(\d+)',dset)
#                 if REG: 
#                     infos['rl']=int(REG.groups()[0])
#                 REG = re.search('c=(\d+)',dset)
#                 if REG: 
#                     infos['c']=int(REG.groups()[0])

#                 subgrid = f[dset]
#                 try:
#                     delta = subgrid.attrs['delta']
#                     origin = subgrid.attrs['origin']
#                     size = subgrid.shape
#                     dim = len(size)
#                     coord = ['x', 'y', 'z']
#                     for i in range(dim) :
#                         infos[coord[i]] = np.arange(0,size[(dim-1)-i])*delta[i]+origin[i]
#                 except:
#                     print(dset)
#                 infos['data'] = np.array(subgrid) 
#                 p.append(infos)

#     return p


# def hdf5_2d(X, Y, data, title=None, colormap='RdBu'):
#     """
#     Create a pseudocolor plot

#     .. note::

#         The dimensions of X and Y should be one greater than those of C. Alternatively, X, Y and C may have equal dimensions, in which case the last row and column of C will be ignored.

#     :param array X: A 1-D array. They will be expanded as needed into the appropriate 2-D arrays, making a rectangular grid.
#     :param array Y: A 1-D array. They will be expanded as needed into the appropriate 2-D arrays, making a rectangular grid.
#     :param array data: A scalar 2-D array. The values will be color-mapped.
#     :param str title: Set figure title.
#     :param str colormap: A Colormap name. The colormap maps the C values to colors. 
#     """
#     size = data.shape
#     fig, ax = plt.subplots()
#     tmpX, tmpY = np.meshgrid(X, Y)
#     im = plt.pcolormesh(tmpX, tmpY, data, cmap=colormap)
#     plt.colorbar(im)
#     ax.set_title(title)
#     plt.show()