"""
`outputdir` get all data file in a simulation directory.
"""

import os
import re


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

def rm_output_active(files):
    """
    Args:
        files (list): A list of file in absolute path

    Return:
        A list of file, not contain output-\d\d\d\d-active
    """
    active_pat = re.compile("\S*/output-(\d\d\d\d)-active/\S*")
    return [f for f in files if not active_pat.match(f)]

def filter_file(files, file_style):
    """
    Choose some file end with file_style

    * par: parameter file
    * scalar file
    * ASCII file
    * HDF5 file
    * checkpoints
    """
    re_pat = {
        "par": "\S*\.par",
        "scalar": "\S*\.(minimum|maximum|norm1|norm2|norm_inf|average)?\.asc(\.(gz|bz2))?$",
        "asc": "\S*\.[xyz]*\.asc(\.(gz|bz2))?$",
        "hdf5": "\S*\.[xyz]*\.h5(\.(gz|bz2))?$",
        "checkpoints" : "\S*/checkpoints\S*",
        "debug": "\S*NaNmask\.\S*\.h5"
    }
    return [f for f in files if re.compile(re_pat[file_style]).match(f)]