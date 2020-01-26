from scipy.interpolate import RegularGridInterpolator

class Interpolator:
    """
    Class representing a grid function defined on a Carpet grid hierarchy.

    :param str interp : interpolation type. 

        * None: no interpolation. The gridarray can only be evaluated at grid points.
        * nearest: nearest neighbour interpolation
        * linear: linear interpolation

    :param values: tuple of nd.array of float, with shapes (m1, ), ..., (mn, ). The points defining the regular grid in n dimensions.
    :param data: array_like, shape (m1, ..., mn, ...). The data on the regular grid in n dimensions.
    """
    def __init__(self, points, values):
        self.grid = tuple([np.asarray(p) for p in points])
        self.values = values