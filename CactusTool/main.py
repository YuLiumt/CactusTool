"""
This is the main modular. Everything start from it.
"""

from .funcs import *
from .parfile import ParFile
from .Carpet import CarpetIOScalar
from .Carpet import CarpetIOHDF5
from .Carpet import CarpetIOASCII
from .debug import DeBug
import os

class Simulation:
    """
    Basis class of CactusTool, anything start from it. Please use it attributes.
    """
    def __init__(self, simpath):
        """
        Args:
            simpath (str): absolute path to simulation directory.
        """
        assert os.path.exists(simpath), "{} doesn't exist in your local computer.".format(simpath)
        
        self.basedir, self.simname = os.path.split(simpath)

        files_tem = fetch_all_file(simpath)

        self.allfiles = rm_output_active(files_tem) # Exclude file in output-0000-active directory.
        self.parfiles = filter_file(self.allfiles, "par")
        self.scafiles = filter_file(self.allfiles, "scalar")
        self.ascfiles = filter_file(self.allfiles, "asc")
        self.h5files  = filter_file(self.allfiles, "hdf5")
        self.debugfiles = filter_file(self.allfiles, "debug")
 
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
            return CarpetIOScalar(self.scafiles)
        else:
            raise Exception("No Scalar variable in {}:".format(self.simname))

    @property
    def H5(self):
        if bool(self.h5files):
            return CarpetIOHDF5(self.h5files)
        else:
            raise Exception("No H5 variable in {}:".format(self.simname))

    @property
    def ASCII(self):
        if bool(self.ascfiles):
            return CarpetIOASCII(self.ascfiles)
        else:
            raise Exception("No ASCII variable in {}:".format(self.simname))        

    @property
    def Debug(self):
        if bool(self.debugfiles):
            return DeBug(self)
        else:
            raise Exception("No NaNCheck in {}:".format(self.simname))