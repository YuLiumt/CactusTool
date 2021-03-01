from ..funcs.file import columns_asc
from ..funcs.log import logger
import numpy as np
import os
import re


class CarpetIOScalar:

    def __init__(self, files):
        """
        Thorn CarpetIOScalar provides I/O methods for outputting scalar values in ASCII format into files. This modular processes output of CarpetIOScalar.

        :param list files: Absolute path to the data file.
        """
        self.files = files
        self.fname = set()
        pat = re.compile('^([a-zA-Z0-9_-]+)\.(minimum|maximum|norm1|norm2|average)?\.asc$')
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
        fileList = [f for f in self.files if key in f]
        return Scalar(fileList)

    def __contains__(self, key):
        return key in self.fname


class Scalar:
    """
    For scalar, we don't need consider grid structure.

    :param list files: Absolute path to the data file.
    """

    def __init__(self, files):
        self.files = files
        self.columns = columns_asc(files[0])
        for file in files:
            assert self.columns == columns_asc(file), "Check why the columns are different."

    
    @property
    def vars(self):
        itime = self.columns.index('time')
        if itime == 1:
            return self.columns[2:]
        elif itime == 8:
            return self.columns[12:]
        else:
            raise Exception("File: {} Header fail to identify.".format(self.files[0]))

    def dsets(self, var='data'):
        p = np.array([[], []])

        itime = self.columns.index('time')
        ivar = self.columns.index(var)

        for f in self.files:
            tem = np.loadtxt(f, usecols=(itime, ivar), comments="#", unpack=True)
            p = np.append(p, tem, axis=1)

        # the sorted unique array
        return np.unique(p, axis=1)