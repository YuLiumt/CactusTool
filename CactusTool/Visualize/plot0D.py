"""
This module provides various functions to plot time series.
"""
import matplotlib.pyplot as plt
from ..utils.units import UnitConversion
import numpy as np

def TimeSeries(ax, dsets, scale=None, unit=('1', '1'), label=None, **kwargs):
    """
    Plot time series plot.

    :param ax: Which axes to use for the plot.
    :param dsets: 2-D array
    """

    t = dsets[0] * UnitConversion(unit[0])
    y = dsets[1] * UnitConversion(unit[1])
    if scale is None:
        im = ax.plot(t, y, **kwargs)
    elif scale == 'log':
        im = ax.semilogy(t, y, **kwargs)
    ax.tick_params(axis='both', direction='in')
    ax.set_xlim(t.min(), t.max())
    if label:
        ax.set_xlabel(label[0])
        ax.set_ylabel(label[1])

    return im