from .check import *
import h5py
import gzip
import bz2
import re
import os


def read(file):
    """
    Carpet output file type is completely different. This function provides different way to open it.

    :param str file: file in absolute path

    >>> with read(file) as f:
    """
    assert os.path.exists(file), "{} doesn't exist in your local computer".format(file)

    if file.endswith('.bz2'):
        f = bz2.BZ2File(file)
    elif file.endswith('.gz'):
        f = gzip.GzipFile(file)
    elif file.endswith('.h5'):
        f = h5py.File(file, 'r')
    elif file.endswith(('.par', '.asc', '.txt')):
        f = open(file)
    else:
        raise RuntimeError("CarpetDataset can't handle this *{} type of file".format(os.path.splitext(file)[1]))
    return f

def fetch_all_file(path):
    """
    Fetch all important file in path, especially `.par`, `.asc`, and `.h5` file.

    :param str path: absolute path
    :return: file list with absolute path.
    """
    assert os.path.exists(path), "{} doesn't exist in your local computer.".format(path)

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
    output-\d\d\d\d-active/ is a copy from output-\d\d\d\d/, so we remove these file in output-\d\d\d\d-active/ to avoid duplicate.

    :param list path: A list of file in absolute path.
    :return: file list withnot the one in output-\d\d\d\d-active
    """
    # Avoid deal with empty files
    is_empty(files)

    files = ensure_list(files)
    active_pat = re.compile("\S*/output-(\d\d\d\d)-active/\S*")
    return [f for f in files if not active_pat.match(f)]

def filter_file(files, file_style):
    """
    Choose the file end with specified file_style

    :param list path: A list of file in absolute path.
    :param str file_style: There are few file_style you can choose:

        * par: parameter file
        * scalar file
        * ASCII file
        * HDF5 file
        * checkpoints
    """
    # Avoid deal with empty files
    is_empty(files)

    files = ensure_list(files)
    re_pat = {
        "par": "\S*\.par",
        "scalar": "\S*\.(minimum|maximum|norm1|norm2|norm_inf|average)?\.asc(\.(gz|bz2))?$",
        "asc": "\S*\.[xyz]*\.asc(\.(gz|bz2))?$",
        "hdf5": "\S*\.[xyz]*\.h5(\.(gz|bz2))?$",
        "checkpoints" : "\S*/checkpoints\S*",
        "debug": "\S*NaNmask\.\S*\.h5"
    }
    return [f for f in files if re.compile(re_pat[file_style]).match(f)]

def columns_asc(file):
    """
    Fetch ASCII file header information.

    :param str file: file in absolute path
    :return: The columns in given file.
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
