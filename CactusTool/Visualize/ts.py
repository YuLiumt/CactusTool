"""
This module provides various functions to plot time series.
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


class TimeSeries:
    """
    Plot time series plot.

    :param ax: Which axes to use for the plot.
    :param ds: Dataset
    """
    def __init__(self, ax, ds):
        self.ax = ax
        if isinstance(ds, pd.DataFrame):
            assert 'time' in ds.columns, "time not in dataset columns"
            self.ds = ds
        
    def plot(self, var, **kwargs):
        im = self.ax.plot(self.ds['time'], self.ds[var], **kwargs)
        self.ax.set_xlabel('time (M)')
        return im
    
    # def annotate(self):
    #     ax.annotate('data = (%.1f, %.1f)'%(xdata, ydata),
        #             (xdata, ydata), xytext=(-2*offset, offset), textcoords='offset points',
        #             bbox=bbox, arrowprops=arrowprops)