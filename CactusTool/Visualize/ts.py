"""
This module provides various functions to plot time series.
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from ..funcs import *


class TimeSeries:
    """
    Plot time series plot.

    :param ax: Which axes to use for the plot.
    :param ds: Dataset
    """
    def __init__(self, ax, t, y):
        self.ax = ax
        self.t = t
        self.y = y
        
    def plot(self, **kwargs):
        """
        line plot.

        :param str fields: The name of the field(s) to be plotted.
        """
        im = self.ax.plot(self.t, self.y, **kwargs)
        self.ax.set_xlabel('time (M)')
        return im
