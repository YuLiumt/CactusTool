"""
This modular processes dataset produced by `Carpet <https://carpetcode.org>`_, which is an adaptive mesh refinement and multi-patch driver for the Cactus.
"""


class CarpetDataset:

    def __init__(self, file):
        self.file = file
        self.filename = None
