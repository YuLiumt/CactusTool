"""
A rectangular grid, supporting 
common arithmetic operations. The standard mathematical unary
functions, e.g. sin(), are implemented as member functions.
Supports interpolation and resampling. The objects can also be used 
as a function, mapping coordinate to interpolated data. Supports 
numercial differentiation.
:ivar data: The actual data (numpy array).
"""

class SlicePlot:
    def __init__(self, ax, ds):
        self.ax = ax
        self.ds = ds

    # def extent():
    #     if extent[0] is not None:
    #         self.ax.set_xlim(xmin=extent[0])
    #     if extent[1] is not None:
    #         self.ax.set_xlim(xmax=extent[1])
    #     if extent[2] is not None:
    #         self.ax.set_ylim(ymin=extent[2])
    #     if extent[3] is not None:
    #         self.ax.set_ylim(ymax=extent[3])