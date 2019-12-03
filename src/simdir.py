from CactusTool.parfile import load_parfile
from CactusTool.CarpetScalars import ScalarsDir
import os

def lazy_property(fn):
    attr_name = '_lazy_' + fn.__name__
    @property
    def _lazy_property(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fn(self))
        return getattr(self, attr_name)
    return _lazy_property

class SimDir:
    """
    This class represents a CACTUS simulation directory.
    """
    def __init__(self, sim, BaseDir='/Users/liuyu/simulations'):
        """
        param:
            path: Path to simulation directory.
        """
        self.path = os.path.join(BaseDir, sim)
        if not os.path.isdir(self.path):
            raise RuntimeError("Folder does not exist: %s" % path)

        self._scan_folders()

    def _scan_folders(self):
        excludes = set(['SIMFACTORY'])
        
        self.dirs     = []
        self.parfiles = []
        self.logfiles = []
        self.errfiles = []
        self.allfiles = []

        def listdir(path):
            l = [os.path.join(path, p) for p in os.listdir(path)]
            return l

        def filter_ext(files, ext):
            return [f for f in files if os.path.splitext(f)[1] == ext]

        def walk_rec(path, level=0):
            self.dirs.append(path)
            if level >= max_depth:
                return
            a = listdir(path)
            f = filter(os.path.isfile, a)
            d = filter(os.path.isdir, a)
            self.allfiles += f
            for p in d:
                if os.path.isdir(p) and (os.path.basename(p) not in excludes):
                    walk_rec(p, level+1)

        max_depth = 4
        walk_rec(self.path)

        self.logfiles = filter_ext(self.allfiles, '.out')
        self.errfiles = filter_ext(self.allfiles, '.err')
        self.parfiles = filter_ext(self.allfiles, '.par')
        self.parfiles.sort(key=os.path.getmtime)
        self.logfiles.sort(key=os.path.getmtime)
        self.errfiles.sort(key=os.path.getmtime)

        if bool(self.parfiles):
            self.params = load_parfile(self.parfiles[0])
        else:
            print("Do not find any parfile in ", self.path)

    @lazy_property
    def Scalars(self):
        return ScalarsDir(self)