"""
Cactus dataset main produced by `Carpet <https://carpetcode.org>`_, which is an adaptive mesh refinement and multi-patch driver for the Cactus. This modular processes output of CarpetIOHDF5.
"""

from ..funcs import *
from multiprocessing import Pool
import numpy as np
import os
import re

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

    @property
    def header_info(self):
        p = {}
        with Pool() as pool:
            result = pool.map(header_h5, self.files)
        for item in result:
            p.update(item)
        return p

    def __str__(self):
        fileList = set(os.path.basename(file) for file in self.files)
        return "\n".join(fileList)