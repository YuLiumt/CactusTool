"""
This is the main modular. Everything start from it.
"""

from .funcs.log import logger
from .funcs.file import is_simfactory
# fetch_all_file, rm_dir_in_filelist, filter_file, 
# from .Parameter import ParFile
from .Carpet.Scalar import CarpetIOScalar
from .Carpet.GF import CarpetGF
# from .Carpet import CarpetIOASCII
# from .Debug import NaNCheck
# from loguru import logger
import os
import re
import glob


class load:

    def __init__(self, simname, basedir=None):
        """
        Load Cactus simulation name.

        :param str simname: simulation name.
        :param str basedir: Basis directory. The default is '~/simulations/'.
        """
        self.simname = simname
        # Make sure basis directory exists.
        if basedir is None:
            self.SourceDir = os.path.expanduser('~/simulations/')
            assert os.path.exists(self.SourceDir), "The default basis directory ~/simulations/ not exists.'"
        else:
            if not basedir.endswith('/'):
                basedir += '/'
            self.SourceDir = basedir
            assert os.path.exists(self.SourceDir), "Basis directory '{}' not exists.".format(self.SourceDir)
        # Make sure simulation directory exists.
        self.simpath = self.SourceDir + self.simname
        assert os.path.exists(self.simpath), "simulation name '{}' not in your '{}'.".format(self.simname, self.SourceDir)
        if not self.simpath.endswith('/'):
            self.simpath += '/'
        # The directory structure of SimFactory is different from the traditional one.
        self._simfactory = is_simfactory(self.simpath)

        """
        # Fetch all file under directory of simpath
        self._allfiles = fetch_all_file(self.simpath)

        # CactusTool is currently unable to deal with the following folder.
        if self._simfactory:
            self._allfiles = rm_dir_in_filelist(self._allfiles, 'SIMFACTORY')
            self._allfiles = rm_dir_in_filelist(self._allfiles, 'output-(\d\d\d\d)-active')
            self._allfiles = rm_dir_in_filelist(self._allfiles, 'cactus-source')
        self._allfiles = rm_dir_in_filelist(self._allfiles, 'checkpoints')
        """

    # def Parfile(self, file=None):
    #     """
    #     Load parameter file if it exist. You can change the default file by use ‘Parfile(<parameter file>)’.

    #     :param str file: parameter file in absolute path.
    #     """
    #     if file:
    #         assert os.path.exists(file), "parameter file '{}' not exists. Make sure it‘s an absolute path.".format(file)
    #         assert file.endswith('.par'), "parameter file '{}' should end with '.par'.".format(file)
    #         self.parfile = file
    #     else:
    #         parfiles = filter_file(self._allfiles, "parfile")

    #         # In some cases, it may contain more than one parfile.
    #         if len(parfiles) == 1:
    #             self.parfile = parfiles[0]
    #         else:
    #             # Guess parfile you want to load
    #             if self._simfactory:
    #                 self.parfile = self.simpath + '/output-0000/' + self.simname + '.par'
    #             else:
    #                 self.parfile = self.simpath + '/' + self.simname + '.par'
    #             assert self.parfile in parfiles, "Make sure `IO::out_dir = $parfile` in your parfile, or you can input the one you want to load."

    #         logger.info("Use the default parameter file '{}'.", self.parfile)

    #     return ParFile(self.parfile)

    def Scalar(self, type='maximum'):
        assert type in ['', 'minimum', 'maximum', 'norm1', 'norm2', 'average']
        if self._simfactory: 
            raise Exception("CactusTool currently cannot handle simfactory!")
        else:
            fileList = glob.glob(self.simpath + '*.{}.asc'.format(type))
        assert bool(fileList), "{} don't have '{}' operation on scalar".format(self.simname, type)
        return CarpetIOScalar(fileList)  
        
    def GF(self, dim='xy', format='h5'):
        assert dim in ['x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz']
        assert format in ['asc', 'h5']
        if self._simfactory: 
            raise Exception("CactusTool currently cannot handle simfactory!")
        else:
            fileList = glob.glob(self.simpath + '*.{}.{}'.format(dim, format))
        assert bool(fileList), "{} don't have '{}' dim in '.{}' format".format(self.simname, dim, format)
        return CarpetGF(fileList, dim, format)

    # @property
    # def ASCII(self):
    #     self.ascfiles = filter_file(self._allfiles, "ascii")
    #     if bool(self.ascfiles):
    #         return CarpetIOASCII(self.ascfiles)
    #     else:
    #         raise Exception("No ASCII variable in {}:".format(self.simname))        

    def Analysis(self, Thorn):
        """
        Analysis thorn's output.

        :param str thorn: thorn name
        """
        from .Analysis import ThornFile, hmns, multipole, volumeintegrals_grmhd

        thorn = Thorn.lower()
        assert thorn in ThornFile.keys(), "CactusTool currently not support {}".format(Thorn)
        if self._simfactory: 
            raise Exception("CactusTool currently cannot handle simfactory!")
        else:
            fileList = glob.glob(self.simpath + ThornFile[thorn])
        assert bool(fileList), "There are no data files about {}".format(Thorn)
        return locals()[thorn](fileList)

    # @property
    # def NaN(self):
    #     self.debug  = filter_file(self.allfiles, "debug")
    #     if bool(self.debugfiles):
    #         return NaNCheck(self)
    #     else:
    #         raise Exception("No NaNCheck in {}:".format(self.simname))