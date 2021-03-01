from .funcs.log import logger
from .funcs.file import is_simfactory
from .funcs.check import ensure_list
import os
# import re
import glob


class load:
    """
    This is the main modular. Everything start from it.
    """

    def __init__(self, simname, basedir='~/simulations/', output=-1):
        """
        Load a Einstein Toolkit simulation.

        :param str simname: simulationname
        :param str basedir: basedir
        :param output: Specify the desired output segments.
        """
        if simname.endswith(os.path.sep):
            simname = simname[:-1]
        self.simname = simname
        basedir = os.path.expanduser(basedir)
        self.simpath = os.path.join(basedir, simname)
        # Make sure simulation directory exists.
        assert os.path.exists(self.simpath), "simulation name '{}' does not exist at path '{}'.".format(self.simname, basedir)
        # The directory structure of SimFactory is different from the traditional one.
        self._simfactory = is_simfactory(self.simpath)
        # output number
        if self._simfactory:
            if output == -1:
                self.output = [i for i in os.listdir(self.simpath) if i[-4:].isdigit()]
            else:
                self.output = ensure_list(output)


    def Scalar(self, ftype='maximum'):
        """
        CarpetIOScalar output

        :param str ftype: reduction operation
        """
        from .Carpet.Scalar import CarpetIOScalar

        assert ftype in ['', 'minimum', 'maximum', 'norm1', 'norm2', 'average']

        if self._simfactory:
            # Data file is stored in multiple folders for SimFactory.
            fileList = []
            for n in self.output:
                path = os.path.join(self.simpath, n, self.simname)
                fileList += glob.glob(path+os.path.sep+'*.{}.asc'.format(ftype))
        else:
            fileList = glob.glob(self.simpath+os.path.sep+'*.{}.asc'.format(ftype))
        return CarpetIOScalar(fileList)


    def GF(self, dim='xy', ftype='h5'):
        """
        CarpetIOHDF5 or CarpetIOASCII output

        :param str dim: dimension
        :param str ftype: endwith
        """
        from .Carpet.GF import CarpetGF
        assert dim in ['x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz']
        assert ftype in ['asc', 'h5']
        if self._simfactory: 
            fileList = glob.glob(self.simpath+'output-????/'+self.simname+'/*.{}.{}'.format(dim, ftype))
        else:
            fileList = glob.glob(self.simpath + '*.{}.{}'.format(dim, ftype))
        assert bool(fileList), "{} don't have '{}' dim in '.{}' ftype".format(self.simname, dim, ftype)
        return CarpetGF(fileList, dim, ftype)


    def Parfile(self, file=None):
        """
        Load parameter file if it exist. You can change the default file by use ‘Parfile(<parameter file>)’.

        :param str file: parameter file in absolute path.
        """
        from .Parameter.parfile import ParFile

        if file:
            assert os.path.exists(file), "parameter file '{}' not exists. Make sure it‘s an absolute path.".format(file)
            assert file.endswith('.par'), "parameter file '{}' should end with '.par'.".format(file)
            self.parfile = file
        else:
            if self._simfactory: 
                raise Exception("CactusTool currently cannot handle simfactory!")
            else:
                fileList = glob.glob(self.simpath + '*.par')
                if len(fileList) == 1:
                    self.parfile = fileList[0]
                else:
                    # Guess parfile you want to load
                    self.parfile = self.simpath + self.simname + '.par'
                    assert self.parfile in fileList, "Make sure `IO::out_dir = $parfile` in your parfile, or you can input the one you want to load."

            logger.info("Use the default parameter file '{}'.", self.parfile)
        return ParFile(self.parfile)


    def Analysis(self, Thorn, fname=None):
        """
        Analysis thorn's output. 

        The current support is as follows:
        [multipole, volumeintegrals_grmhd, puncturetracker]

        :param str thorn: thorn name
        :param str fname: file name
        """
        from . import Analysis

        thorn = Thorn.lower()

        ThornFile = {
            'hmns': 'HMNS_*.asc',
            'multipole': 'mp_*',
            'volumeintegrals_grmhd': 'volume_integrals-GRMHD.asc',
            'twopunctures': 'TwoPunctures.bbh',
            'puncturetracker': 'puncturetracker-pt_loc..asc',
        }
        if fname is None:
            assert thorn in ThornFile.keys(), "Use fname param"
            fname = ThornFile[thorn]
        
        if self._simfactory:
            # Data file is stored in multiple folders for SimFactory.
            fileList = []
            for n in self.output:
                path = os.path.join(self.simpath, n, self.simname)
                fileList += glob.glob(path+os.path.sep+fname)
        else:
            fileList = glob.glob(self.simpath+os.path.sep+fname)
        assert bool(fileList), "There are no data files about {}".format(Thorn)

        try:
            return getattr(Analysis, thorn)(fileList)
        except AttributeError:
            print("CactusTool currently not support %s thorn" % (Thorn))
