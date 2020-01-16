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