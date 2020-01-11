"""
This modular provided many useful functions needed by other modular.
"""

import numpy as np
import json
import h5py
import gzip
import bz2
import re
import os


###################    file    ###################

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


###################    python type    ###################

def is_empty(obj):
    """
    Make sure `obj` isn't empty.
    """
    if not obj:
        print("Please check why {} is None".format(obj))

def ensure_list(obj):
    """
    This function ensures that `obj` is a list. Typically used to convert a string to a list.
    """
    if obj is None:
        return [obj]
    if not isinstance(obj, list):
        return [obj]
    return obj

def ensure_numpy_array(obj):
    """
    This function ensures that `obj` is a numpy array. Typically used to convert scalar, list or tuple argument passed to functions using Cython.
    """
    if isinstance(obj, np.ndarray):
        if obj.shape == ():
            return np.array([obj])
        # We cast to ndarray to catch ndarray subclasses
        return np.array(obj)
    elif isinstance(obj, (list, tuple)):
        return np.asarray(obj)
    else:
        return np.asarray([obj])

###################    json    ###################

def Format(dicts):
    return json.dumps(dicts, sort_keys=True, indent=4)

def subkey_have_value(dicts, subkey, value):
    p = {}
    for k, v in dicts.items():
        assert subkey in v, "dict don't have subkey: %s" % (subkey)
        if v[subkey] == value:
            p[k] = dicts[k]
    return p

def subkey_contain_element(dicts, subkey):
    p = set()
    for k, v in dicts.items():
        assert subkey in v, "dict don't have subkey: %s" % (subkey)
        try:
            p.add(v[subkey])
        except:
            raise RuntimeError("the value of subkey is not element" % subkey)
    return p



def add_two_key(par, key_a, key_b, val):
    if key_a in par:
        par[key_a].update({key_b: val})
    else:
        par.update({key_a: {key_b: val}})

def add_three_key(par, key_a, key_b, key_c, val):
    if key_a in par:
        if key_b in par:
            par[key_a][key_b].update({key_c: val})
        else:
            par[key_a].update({key_b: {key_c: val}})
    else:
        par.update({key_a: {key_b: {key_c: val}}})