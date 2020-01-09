"""
This module provides functions to plot different type of figure.

* color_bar: :py:func:`color_bar`
* 2-D scalar plot: :py:func:`imshow`, :py:func:`pcolormesh`, :py:func:`contour`, :py:func:`surface`
* 2-D vector plot: :py:func:`vectors`
* grid struct: :py:func:`grid_struct`
"""

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import numpy as np


def color_bar(image, **kwargs):
    """
    Adds a colorbar with given axes. 

    :param im: Which image object to use for the colorbar.
    :param kwargs: Unknown keyword arguments are passed to :py:func:`plt.colorbar()`.
    """
    plt.colorbar(image, **kwargs)


def imshow(ax, Z, axis_ticks=False, **kwargs):
    """
    Plot an array as an image. Loss coordinate information. This function return image object, you can add colorbar use :py:func:`color_bar`.

    If X and Y are each equidistant, imshow can be more convenient.

    :param ax: Which axes to use for the plot.
    :param array Z: 2D data to plot.
    :param bool axis_ticks: Whether turn on axis lines and labels.
    :param kwargs: Unknown keyword arguments are passed to :py:func:`ax.imshow()`
    :return: image object

    Supported interpolation are ['none', 'nearest', 'bilinear', 'bicubic', 'spline16', 'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric', 'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos'].

    >>> im = imshow(ax, Z, axis_ticks=True)
    turn on axis lines and labels.
    >>> im = imshow(ax, Z, cmap='RdBu')
    colormap name used to map scalar data to colors.
    >>> im = imshow(ax, Z, alpha=0.5)
    The alpha blending value, between 0 (transparent) and 1 (opaque).
    >>> im = imshow(ax, Z, interpolation='nearest')
    Interpolation method.
    """
    im = ax.imshow(Z, **kwargs)
    if not axis_ticks:
        ax.set_xticks([])
        ax.set_yticks([])
    return im

def pcolormesh(ax, X, Y, Z, **kwargs):
    """
    Create a pseudocolor 2-D plot.

    :param ax: Which axes to use for the plot.
    :param array X, Y: The coordinates of X and Y are bounds.
    :param array Z: 2D data to plot. 
    :param kwargs: Unknown keyword arguments are passed to :py:func:`ax.pcolormesh()`.
    :return: image object

    .. note::

        If X, Y and Z have equal dimensions, the last value from the z array may be removed.

    >>> im = pcolormesh(ax, x, y, Z)
    Create a pseudocolor 2-D plot.
    >>> Zm = np.ma.masked_where(np.abs(Z) < 0.0001, Z)
    >>> im = pcolormesh(ax, x, y, Zm, shading='gouraud')
    """
    if X.shape == Z.shape and Y.shape == Z.shape:
        print("Warning: The last value from the z array was removed.") 
    im = ax.pcolormesh(X, Y, Z, **kwargs)
    return im

def contour(ax, X, Y, Z, fill=True, contour_levels=10, line_levels=5, line_color='k', **kwargs):
    """
    Draw contour lines or filled contours.

    :param ax: Which axes to use for the plot.
    :param array X, Y: The coordinates of the values in Z. X and Y must both be 2-D with the same shape as Z (e.g. created via :py:func:`numpy.meshgrid`)
    :param array Z: 2D data to plot. 
    :param bool fill: draw contour lines or filled contours
    :param contour_levels: if fill is True; contour_levels determines the number of the filled contours.
    :param line_levels: if fill is True; line_levels determines the number of the lines.
    :param line_color: if fill is True; line_color determines the color of the lines.
    :param kwargs: Unknown keyword arguments are passed to :py:func:`ax.contour()` or :py:func:`ax.contourf()`.
    :return: image object

    >>> im = imshow(ax, X, Y, Z)
    Only draw contour lines.
    >>> im = contour(ax, x, y, Z, fill=True)
    Draw contour lines and filled contours,
    >>> im = contour(ax, x, y, Z, fill=True, contour_levels=10, line_levels=5)
    Draw contour lines and filled contours, filled levels determined by contour_levels and line number determined by line_levels
    >>> im = contour(ax, x, y, Z, fill=True, line_color='r')
    line color determined by line_color
    >>> im = contour(ax, x, y, Z, fill=True, cmap=plt.cm.bone)
    filled color determined by cmap
    >>> im = contour(ax, X, Y, Z, fill=True, alpha=0.5)
    The alpha blending value, between 0 (transparent) and 1 (opaque).
    >>> Zm = np.ma.masked_where(np.abs(Z) < 0.000001, Z)
    >>> im = contour(ax, x, y, Zm, fill=False)
    control the masked region.
    >>> from matplotlib import ticker
    >>> im = contour(ax, x, y, Z, locator=ticker.LogLocator(), fill=True)
    log locator tells contourf to use a log scale
    """
    assert X.shape == Z.shape and Y.shape == Z.shape, "The coordinates of X and Y must have the same shape as Z"
    if fill:
        line_factor = int(contour_levels)//int(line_levels)
        im = ax.contourf(X, Y, Z, levels=contour_levels, **kwargs)
        ax.contour(im, levels=im.levels[::line_factor], colors=line_color)
    else:
        im = ax.contour(X, Y, Z, locator=plt.LogLocator(), **kwargs)
        fmt = ticker.LogFormatterMathtext()
        ax.clabel(im, im.levels, fmt=fmt)
    return im

def surface(ax, X, Y, Z, mesh=False, **kwargs):
    """
    Create a surface plot.

    :param ax: Which axes to use for the plot. ax is the 3D projection.
    :param array X, Y: The coordinates of the values in Z. X and Y must both be 2-D with the same shape as Z (e.g. created via :py:func:`numpy.meshgrid`)
    :param array Z: 2D data to plot.  
    :param bool mesh: A mesh plot
    :param kwargs: Unknown keyword arguments are passed to :py:func:`ax.contour()` or :py:func:`ax.contourf()`.
    :return: image object

    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111, projection='3d')
    >>> im = surface(ax, x, y, Z)
    Create a surface plot.
    >>> im = surface(ax, x, y, Z, cmap=plt.cm.YlGnBu_r)
    A colormap for the surface patches.
    >>> im = surface(ax, x, y, Z, mesh=True)
    A mesh plot.
    """
    if mesh:
        im = ax.plot_wireframe(X, Y, Z, **kwargs)
    else:
        im = ax.plot_surface(X, Y, Z, **kwargs)
    return im

def vectors(ax, X, Y, U, V, W=None, **kwargs):
    """
    Plot vectors in 2-D plot. The color of the arrows is its mod. :py:func:`ax.quiverkey()` auto-scales the length of the arrows to a reasonable size. 

    :param ax: Which axes to use for the plot.
    :param array X, Y: X, Y define the arrow locations.
    :param array U, V, W: 3-D vector have component (U, V, W). Default is 2-D vector (U,V)
    :param kwargs: Unknown keyword arguments are passed to :py:func:`ax.quiverkey()`.
    :return: image object

    >>> im = vectors(ax, x, y, U, V)
    Plot a 2D field of vectors.
    >>> im = vectors(ax, x, y, U, V, W)
    Plot a 3D field of vectors in 2D plot.
    """
    M = np.hypot(U, V, W)
    im = ax.quiver(X, Y, U, V, M, **kwargs)
    return im

def grid_struct():
    """
    Plots grid structure.
    """
    im = None
    return im