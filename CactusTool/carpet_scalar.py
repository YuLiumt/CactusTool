"""
The :py:mod:`~.cactus_scalars` module provides functions to load 
timeseries in Cactus formats and a class :py:class:`ScalarsDir` for easy 
access to all timeseries in a Cactus simulation directory. This module 
is normally not used directly, but from the :py:mod:`~.simdir` module. 
The data loaded by this module is represented as 
:py:class:`~.TimeSeries` objects.
"""

import numpy as np
import json
import os
import re

def AllVars(files):
    """
    Args:
        files (list): A list of file in absolute path
    
    Return:
        A disc contain Useful infomation
    """
    for file in files:
        print(os.path.basename(file))
        print(column_header(file))

def column_header(file):
    return None

def variables(files):
    """
    Args:
        files (list): A list of file in absolute path
    
    Return:
        A dict of variables from filename.
    """
    scalar_pat = re.compile("(\S*)\.(minimum|maximum|norm1|norm2|norm_inf|average)?\.asc(\.(gz|bz2))?$")
    vars = {}
    for file in files:
        name = os.path.basename(file)
        var = scalar_pat.match(name).group(1)
        type = scalar_pat.match(name).group(2)
        if var not in vars.keys():
            vars.update({var: {type: file}})
    return vars

def Column(file):
    """
    Get the columns name in Scalar file.

    Args:
        file (str): open file
    
    Return:
        A list of column name
    """
    with read(file) as f:
        columns=[]
        for line in f.readlines():
            if "# 1:iteration 2:time 3:data" in line:
                columns = columns + line.split()[1:]
            if "# column format:" in line:
                columns = line.split()[3:]
            if "# data columns: " in line:
                del columns[-1]
                columns = columns + line.split()[3:]
                break
    if len(columns) > 0:
        return [name.split(":")[1] for c, name in enumerate(columns)]
    else:
        raise Exception("File: {} Header fail to identify.".format(file))

def Data(file):
    """
    Get the data in Scalar file.

    Args:
        file (str): open file
    
    Return:
        data
    """
    return np.loadtxt(file, comments="#", unpack=True)

def merge_components(filelist):
    """
    Some Variable data locate in different components.

    Args:
        file (list): a list of scalar file with different components.
    
    Return:
        data contain all components.
    """
    data = np.array()
    for i in range(len(filelist)):
        try: 
            tmp = np.loadtxt(filelist[i], comments="#", unpack=True)
            data = np.append(data, tmp, axis=1)
        except:
            print("[ERROR] Unable to load file:", filelist[i])
    return data


class VarReader:
    def __init__(self, var):
        self.var = var
        print(self.var, 'No')


class Scalar:
    """
    Operation in Carpet Scalar file.

    :ivar min:       access to minimum reduction.
    :ivar max:       access to maximum reduction.
    :ivar norm1:     access to norm1 reduction.
    :ivar norm2:     access to norm2 reduction.
    :ivar average:   access to average reduction.

    .. note::
        infnorm is reconstructed from min and max if infnorm itself is not available.
    
    :param sd:  Simulation directory.
    :type sd:   :py:class:`~.SimDir` instance.

    """
    def __init__(self, Sim):
        self.files = Sim.scafiles
        self.allvars = AllVars(self.files)
        # self.variables = variables(self.files)

    def __getitem__(self, var):
        return VarReader(self.variables[var])
    
    def __contains__(self, key):
        return key in self.variables

    def __str__(self):
        return "All Scalar in simulations dortory:\n {}".format(json.dumps(self.variables, sort_keys=True, indent=4))


if __name__ == "__main__":
    from simulation import Sim
    sim = Sim('/Users/liuyu/simulations/TOV_single_vel')
    # sim = Sim('/Users/liuyu/simulations/BH')
    a = sim.Scalars
    print(a.allvars)
    # a['hydro_analysis-hydro_analysis_rho_max_soc']
    # print('hydro_analysis-hydro_analysis_rho_max_loc' in a)