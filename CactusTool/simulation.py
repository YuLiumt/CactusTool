"""
`simulation` provides easy access to CACTUS data files.

A simulation directory is represented by an instance of the  
`Sim` class, which provides access to all supported
data types.
"""

from . import outputdir
import os

class Sim:
    """
    Basis class of CactusTool, anything start from it. Please use it attributes. It first read parameter file if it exist.
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
        self.allfiles = outputdir.rmoutputactive(allfiles) # Exclude file in output-0000-active directory.

        self.parfiles = outputdir.filter_par(self.allfiles)
        self.scafiles = outputdir.filter_scalar(self.allfiles)
        self.ascfiles = outputdir.filter_asc(self.allfiles)
        self.h5files  = outputdir.filter_h5(self.allfiles)

        self.has_parfile = bool(self.parfiles)
        # if bool(self.parfiles):
    #         self.params = load_parfile(self.parfiles[0])
    #     else:
    #         print("Do not find any parfile in ", self.path)

    @property
    def Scalars(self):
        if bool(self.scafiles):
            return carpetscalar.Scalar(self)
        else:
            raise Exception("No Scalar variable in {}:".format(simpath))
