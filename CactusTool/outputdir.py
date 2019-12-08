"""
`outputdir` get all data file in a simulation directory.
"""

import os
import re


def filter_par(files):
    """
    Args:
        files (list): A list of file in absolute path.

    Return:
        A list of parameter file in different iteration.
        
    Usage:
        >>> filter_par(files)
        ['/Users/liuyu/simulations/TOV_single_vel/output-0000/TOV_single_vel.par']
    """
    # par_pat = re.compile("\S*/output-(\d\d\d\d)/[^/]*\.par")
    par_pat = re.compile("\S*\.par")
    return [f for f in files if par_pat.match(f)]

def filter_scalar(files):
    """
    Args: 
        files (list): A list of file in absolute path.

    Return:
        A list of scalar data file.
    """
    scalar_pat = re.compile("\S*\.(minimum|maximum|norm1|norm2|norm_inf|average)?\.asc(\.(gz|bz2))?$")
    # scalar_pat = re.compile("\S*/output-(\d\d\d\d)/\S*\.(minimum|maximum|norm1|norm2|norm_inf|average)?\.asc(\.(gz|bz2))?$")
    return [f for f in files if scalar_pat.match(f)]

def filter_asc(files):
    """
    Args:
        files (list): A list of file in absolute path

    Return:
        A list of ASCII grid data file.
    """
    # asc_pat = re.compile("\S*/output-(\d\d\d\d)/\S*\.[xyz]*\.asc(\.(gz|bz2))?$")
    asc_pat = re.compile("\S*\.[xyz]*\.asc(\.(gz|bz2))?$")
    return [f for f in files if asc_pat.match(f)]

def filter_h5(files):
    """
    Args:
        files (list): A list of file in absolute path

    Return:
        A list of HDF5 grid data file.
    """
    # h5_pat = re.compile("\S*/output-(\d\d\d\d)/\S*\.[xyz]*\.h5(\.(gz|bz2))?$")
    h5_pat = re.compile("\S*\.[xyz]*\.h5(\.(gz|bz2))?$")
    return [f for f in files if h5_pat.match(f)]

def filter_check(files):
    """
    Args:
        files (list): A list of file in absolute path

    Return:
        A list of checkpoints file
    """
    check_pat = re.compile("\S*/checkpoints\S*")
    return [f for f in files if check_pat.match(f)]

def fetch_all_datafile(path):
    """
    Only fetch the absolute path of .par, .asc, and .h5 file.

    Args:
        path (str): absolute path

    Return:
        list: All file in path.
    """
    exclude_dirs = set(['SIMFACTORY']) # Exclude SIMFACTORY directory

    filelist = []
    for root, dirs, files in os.walk(path):
        if os.path.basename(root) not in exclude_dirs:  
            for file in files:
                #TODO Some data file may be end with .bz2 or .gz.
                if file.endswith(('.par', '.asc', '.h5')):
                    filelist.append(os.path.join(root, file))

    return filelist

def rmoutputactive(files):
    """
    Args:
        files (list): A list of file in absolute path

    Return:
        A list of file, not contain output-\d\d\d\d-active
    """
    active_pat = re.compile("\S*/output-(\d\d\d\d)-active/\S*")
    return [f for f in files if not active_pat.match(f)]