import matplotlib.pyplot as plt
import numpy as np


def hdf5_2d(X, Y, data, title=None, colormap='RdBu'):
    """
    Create a pseudocolor plot

    .. note::

        The dimensions of X and Y should be one greater than those of C. Alternatively, X, Y and C may have equal dimensions, in which case the last row and column of C will be ignored.

    :param array X: A 1-D array. They will be expanded as needed into the appropriate 2-D arrays, making a rectangular grid.
    :param array Y: A 1-D array. They will be expanded as needed into the appropriate 2-D arrays, making a rectangular grid.
    :param array data: A scalar 2-D array. The values will be color-mapped.
    :param str title: Set figure title.
    :param str colormap: A Colormap name. The colormap maps the C values to colors. 
    """
    size = data.shape
    fig, ax = plt.subplots()
    tmpX, tmpY = np.meshgrid(X, Y)
    im = plt.pcolormesh(tmpX, tmpY, data, cmap=colormap)
    plt.colorbar(im)
    ax.set_title(title)
    plt.show()