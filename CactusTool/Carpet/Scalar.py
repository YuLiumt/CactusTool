from ..utils.file import columns_asc, dataset_asc
from ..utils.check import ensure_list
from ..Analysis import TimeSeries, DataSeries
from functools import cached_property
import os


class CarpetIOScalar:

    def __init__(self, files, ftype):
        """
        This class provides access to various types of scalar data. It's a dictionary-like object with keys the file name.

        :param list files: Absolute path to the data file.
        :param str ftype: Type of reduction.
            - scalar:  access to grid scalars.
            - minimum: access to minimum reduction.
            - maximum: access to maximum reduction.
            - norm1:   access to norm1 reduction.
            - norm2:   access to norm2 reduction.
            - average: access to average reduction.
            - infnorm: access to inf-norm reduction.
        """
        self.files = ensure_list(files)
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

    def __str__(self):
        l = max(len(i) for i in self.fname)
        ret = f"Reduction type: '{self.endwith[1:-4]}'\n"
        ret += "File Name".ljust(l) + " Variables\n"
        ret += "---------".ljust(l) + " ---------\n"
        for i in self.fname:
            ret += "{} {}\n".format("\033[1m"+i.ljust(l)+"\033[0m", self[i].vars)
        return ret


class Scalar:

    def __init__(self, files):
        """
        Single variable per file or single file per group are supported. 

        :param list files: Absolute path to the data file.
        """
        self.files = files
        self.columns = columns_asc(files[0])
        for file in files[1:]:
            assert self.columns == columns_asc(file), "Variable:\n{}:\n{}\n-----------\n{}:\n{}".format(files[0], self.columns, file, columns_asc(file))
    
    @cached_property
    def vars(self):
        itime = self.columns.index('time')
        if itime == 1:
            return self.columns[2:]
        elif itime == 8:
            return self.columns[12:]
        else:
            raise Exception("File: {} Header fail to identify.".format(self.files[0]))

    def dsets(self, vars='data'):
        itime = self.columns.index('time')
        if isinstance(vars, list):
            ivar = [self.columns.index(v) for v in vars]
            ivar.insert(0, itime)
            dset = dataset_asc(self.files, ivar)
            return DataSeries(dset[0], dset[1:].T, vars)
        else:
            ivar = self.columns.index(vars)
            dset = dataset_asc(self.files, (itime, ivar))
            return TimeSeries(dset[0], dset[1])

    def filesize(self, unit="MB"):
        """
        Return the total size of the given files. Available units B, KB, MB and GB
        """
        units = {"B": 1, "KB": 1024, "MB": 1024 ** 2, "GB": 1024 ** 3}
        assert unit in units, "Invalid unit: expected one of {}".format(list(units.keys()))
        return sum(os.path.getsize(file) for file in self.files) / units[unit]
