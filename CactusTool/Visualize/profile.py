"""
This module provides various functions to plot 1-D profile.
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

class LinePlot:
    """
    Plot 1-D profile.

    :param ax: Which axes to use for the plot.
    :param ds: Dataset
    """
    def __init__(self, ax, ds):
        self.ax = ax
        if isinstance(ds, pd.DataFrame):
            assert 'time' in ds.columns, "time not in dataset columns"
            self.ds = ds

    # def plot(self):