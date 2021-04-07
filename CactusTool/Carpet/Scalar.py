from ..utils.file import columns_asc, dataset_asc
from ..Analysis import DataFrameSeries
# from ..utils.log import logger
import pandas as pd
import os


class CarpetIOScalar:

    def __init__(self, files, ftype):
        """
        Thorn CarpetIOScalar provides I/O methods for outputting scalar values in ASCII format into files. This modular processes output of CarpetIOScalar.

        :param list files: Absolute path to the data file.
        """
        self.files = files
        self.endwith = ".{}.asc".format(ftype)
        self.fname = set(os.path.basename(f).replace(self.endwith, '') for f in self.files)


    def __getitem__(self, key):
        """
        Read the dataset when given specified file name.

        :param str key: file name.
        """
        assert key in self.fname, "File {}{} is not exist".format(key, self.endwith)
        fileList = [f for f in self.files if os.path.basename(f).split(".")[0] == key]
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
            assert self.columns == columns_asc(file), "Variable:\n{}:\n{}\n-----------\n{}:\n{}".format(files[0], self.columns, file, columns_asc(file))
    
    @property
    def vars(self):
        itime = self.columns.index('time')
        if itime == 1:
            return self.columns[2:]
        elif itime == 8:
            return self.columns[12:]
        else:
            raise Exception("File: {} Header fail to identify.".format(self.files[0]))

    def dsets(self, vars='data'):
        if not isinstance(vars, list):
            vars = [vars]
        vars.insert(0, 'time')
        ivar = [self.columns.index(v) for v in vars]
        df = pd.DataFrame(dataset_asc(self.files, ivar).T, columns=vars)
        return DataFrameSeries(df.set_index('time'))

