"""
`simulation` provides easy access to CACTUS data files.

A simulation directory is represented by an instance of the `Sim` class, which provides access to all supported data types.
"""

from . import outputdir
from .parfile import ParFile
from .Carpet.Scalar import CarpetIOScalar
from .Carpet.HDF5 import CarpetIOHDF5
from .debug import DeBug
import os

class Sim:
    """
    Basis class of CactusTool, anything start from it. Please use it attributes.
    """
    def __init__(self, simpath):
        """
        Args:
            simpath (str): absolute path to simulation directory.
        """
        if not os.path.exists(simpath):
            raise Exception("Path: {} do not exists:".format(simpath))
        
        self.basedir, self.simname = os.path.split(simpath)

        allfiles = outputdir.fetch_all_datafile(simpath)

        self.allfiles = outputdir.rm_output_active(allfiles) # Exclude file in output-0000-active directory.
        self.parfiles = outputdir.filter_file(self.allfiles, "par")
        self.scafiles = outputdir.filter_file(self.allfiles, "scalar")
        self.ascfiles = outputdir.filter_file(self.allfiles, "asc")
        self.h5files  = outputdir.filter_file(self.allfiles, "hdf5")
        self.debugfiles = outputdir.filter_file(self.allfiles, "debug")
 
    @property
    def Par(self):
        """
        It first read parameter file if it exist.
        """
        if bool(self.parfiles):
            return ParFile(self.parfiles)
        else:
            print("Do not find any parfile in ", self.path)

    @property
    def Scalar(self):
        if bool(self.scafiles):
            return CarpetIOScalar(self)
        else:
            raise Exception("No Scalar variable in {}:".format(self.simname))

    @property
    def H5(self):
        if bool(self.h5files):
            return CarpetIOHDF5(self)
        else:
            raise Exception("No H5 variable in {}:".format(self.simname))

    @property
    def Debug(self):
        if bool(self.debugfiles):
            return DeBug(self)
        else:
            raise Exception("No NaNCheck in {}:".format(self.simname))