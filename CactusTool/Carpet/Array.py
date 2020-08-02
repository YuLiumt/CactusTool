from ..funcs import read, ensure_list
from ..Analysis import HMNS


class IOArray:

    def __init__(self, files):
        self.files = ensure_list(files)

    def __getitem__(self, thorn):
        """
        Specify the thorn

        :param str thorn: thorn name
        """
        return globals()[thorn](self.files)

