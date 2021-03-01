"""
plot 1D Carpet data at a given iteration.
"""
import matplotlib.pyplot as plt
# from ..funcs.unit import 
import numpy as np

def LinePlot(ax, dsets, fmt=None, reflevel=-1, unit='cactus', **kwargs):
    """
    plot 1D Carpet data at a given iteration.

    :param ax: Which axes to use for the plot.
    :param dict dsets: Dataset
    :param int reflevel: only the refinement level specified in it will be displayed.
    """
    if fmt:
        line = ax.plot([], [], fmt, **kwargs)[0]
    else:
        line = ax.plot([], [], **kwargs)[0]
    for rl in sorted(dsets):
        if reflevel != -1 and rl != reflevel:
            continue
        for c in sorted(dsets[rl]):
            coord = dsets[rl][c]['coord']
            data = dsets[rl][c]['data']
            time = dsets[rl][c]['time']
            # Add the data to the line object
            line.set_xdata(np.append(line.get_xdata(), coord))
            line.set_ydata(np.append(line.get_ydata(), data))
    # Sort points
    indices = np.argsort(line.get_xdata())
    line.set_xdata(line.get_xdata()[indices])
    line.set_ydata(line.get_ydata()[indices])
    # Adjust the axes
    ax.relim()
    ax.autoscale_view()
    if unit == 'cactus':
        ax.set_xlabel('Coord [M]')
    ax.set_title("Time: {}".format(time))
    return line