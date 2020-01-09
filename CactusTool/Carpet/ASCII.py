"""
Cactus dataset main produced by `Carpet <https://carpetcode.org>`_, which is an adaptive mesh refinement and multi-patch driver for the Cactus. This modular processes output of CarpetIOASCII.
"""

class CarpetIOASCII:
    """
    This class can handle all CarpetIOASCII output. Read the dataset from files and return all useful information stored in these files.

    :param list files: A list of file in absolute path.
    """
    def __init__(self, files):
        self.files = ensure_list(files)
        self.x     = Griddim(self.files, 'x')
        self.y     = Griddim(self.files, 'y')
        self.z     = Griddim(self.files, 'z')
        self.xy    = Griddim(self.files, 'xy')
        self.xz    = Griddim(self.files, 'xz')
        self.yz    = Griddim(self.files, 'yz')
        self.xyz   = Griddim(self.files, 'xyz')


class Griddim:
    def __init__(self, files, dim):
        self.dim = dim
        pat_fn = re.compile("\S*\.([xyz]*)\.asc(\.(gz|bz2))?$")
        self.files = [file for file in files if pat_fn.match(file).group(1) == self.dim]