from scipy.interpolate import griddata
import numpy as np
import pandas as pd


class AMRGrid:
    def __init__(self, dsets, dim, var):
        self.dsets = dsets
        self.dim = dim
        self.var = var

    @property
    def hierarchy(self):
        p = {}

        rls = self.dsets['rl'].unique().astype(int).tolist()
        for rl in rls:
            dset = self.dsets[self.dsets.rl == rl]
            c = dset.c.unique().astype(int).tolist()
            p[rl] = c

        return p


    def coords(self, rl=-1, c=-1):
        slice = self.dsets
        if rl != -1:
            slice = slice[slice.rl == rl]
            assert not slice.empty, "dataset is empty at refinement level {}".format(rl)
        if c != -1:
            slice = slice[slice.c == c]
            assert not slice.empty, "dataset is empty at component {}".format(c)
        grid = tuple(np.sort(np.unique(slice[dim].values)) for dim in self.dim)
        return tuple(np.meshgrid(*grid))

    def interpolate(self, coords):
        points = tuple([self.dsets[dim].values for dim in self.dim])
        return griddata(points, self.dsets[self.var].values, coords, method='nearest')






    # def mesh(self):
        # for d in range(self.ndim):
        # L[d] = self.origin[d] + self.delta[d]*L[d]

        # return tuple(L)

    # def coordinate(self, index):
    #     """
    #     Coordinates of the given index array
    #     """
    #     return self.origin + (index - self.iorigin) * self.delta

class UniformGrid:
    def __init__(self, grid):
        self.origin = grid['origin']
        self.delta = grid['delta'] 
        self.reflevel = grid['level']
        self.nghost = grid['cctk_nghostzones'] 
        self.time = grid['time']
        self.data = grid['data']
        self.shape = self.data.shape

    def scale_coords(self, scale):
        """Rescale all coordinates.

        :param scale: Factor to scale by.
        :type scale:  float or 1d numpy array
        """
        self.origin = self.origin * scale
        self.delta  = self.delta * scale

    def shift_coords(self, shift):
        """Shift all coordinates.

        :param shift: vector added to the origin.
        :type shift:  float or 1d numpy array
        """
        self.origin = self.origin + np.asarray(shift)

    def position(self, index):
        return np.array(index) * self.delta + self.origin

    def index(self, position):
        index = (((np.array(position) - self.origin) / self.delta) + 0.5).astype(int)
        return tuple(index)

    @property
    def dv(self):
        """
        :returns: Volume of a grid cell.
        :rtype:   float
        """
        return np.prod(self.delta)

    @property
    def volume(self):
        """
        :returns: Volume of the whole grid.
        :rtype:   float
        """
        return np.prod(self.shape) * self.dv

    # def contains(self, position):
    #     """Test if a coordinate is contained in the grid. The size of the grid cells is taken into account, resulting in a cube larger by dx/2 on each side compared to the one given by x0, x1.

    #     :param pos: Coordinate to test.
    #     :returns:   If pos is contained.
    #     """
    #     if not alltrue( pos > (self.x0() - 0.5 * self.dx()) ):
    #         return False
    #     if not alltrue( pos < (self.x1() + 0.5 * self.dx()) ):
    #         return False
    #     return True

