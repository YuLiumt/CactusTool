from .check import *
import pandas as pd
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
        raise RuntimeError("CactusTool can't handle this *{} type of file".format(os.path.splitext(file)[1]))
    return f

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

def header_h5(file):
    """
    Return a dictionary about header information

    :param str file: file in absolute path.
    :return: a dictionary about header information
    """
    pattern = re.compile(r'([^:]+)::(\S+) it=(\d+) tl=(\d+)( m=(\d+))? rl=(\d+)( c=(\d+))?')
    p = {}
    with read(file) as f:
        for item in list(f):
            dset = pattern.match(item)
            if dset is None:
                continue
            p[item] = {
                'file': file,
                'thorn': dset.group(1), 
                'varname': dset.group(2),
                'iteration': int(dset.group(3)),
                'timelevel': int(dset.group(4)),
                'rl': int(dset.group(7))
            }
            if dset.group(6) != None:
                p[item].update({'ml': int(dset.group(6))})
            if dset.group(9) != None:
                p[item].update({'c': int(dset.group(9))})

    return p

def select_header_h5(header, var=-1, it=-1, rl=-1, c=-1):
    """
    Select specified header in return of :py:func:`header_h5`.

    :param str var: -1: all variable.
    :param int it: -1: all iteration number.
    :param int rl: -1: all refinement level.
    :param int c: -1: all component.
    :return: a list of header
    """
    p = []
    for item in header:
        if var == -1 or header[item]['varname'] == var:
            if it == -1 or header[item]['iteration'] == it:
                if rl == -1 or header[item]['rl'] == rl:
                    if c == -1 or header[item]['c'] == c:
                        p.append(item)

    return p

def dataset_h5(file, slice):
    """
    Return a dictionary about dataset.

    :param str file: file in absolute path.
    :param str slice: header inf
    :return: a dictionary about dataset
    """
    p = dict()
    
    with read(file) as f:
        if slice not in sorted(list(f)):
            raise RuntimeError("%s not in %s" % (slice, file))
        dset = f[slice]
        for item in list(dset.attrs):
            p[item] = dset.attrs.get(item, None)
        p['data'] = np.array(dset)  
    
    return p

def dataframe_h5(file, axis):
    """
    HDF5 to pandas DataFrame
    """
    assert axis in ['x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz'], "Does not include {} axis on {}".format(dim, file)
    pattern = re.compile(r'([^:]+)::(\S+) it=(\d+) tl=(\d+)( m=(\d+))? rl=(\d+)( c=(\d+))?')
    with read(file) as f:
        tem = []
        for item in list(f):
            header = pattern.match(item)
            if header is None:
                continue
            dset = f[item]
            var = header.group(2)
            data = np.array(dset)
            origin = dset.attrs.get('origin', None)
            delta = dset.attrs.get('delta', None)
            size = data.shape
            dim = len(axis)
            coord = tuple(np.arange(0,size[(dim-1)-i])*delta[i]+origin[i] for i in range(dim))
            grid = np.meshgrid(*coord)
            tem_c = pd.DataFrame()
            for i in range(dim):
                tem_c[axis[i]] = grid[i].flatten()
            tem_c[var] = data.flatten()
            tem_c['rl'] = int(header.group(7))
            tem_c['it'] = int(header.group(3))
            tem_c['time'] = dset.attrs.get('time', None)
            if header.group(9) != None:
                tem_c['c'] = int(header.group(9))    
            tem.append(tem_c)

    return pd.concat(tem)

################### Folder ###################
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
    output-\d\d\d\d-active/ is a copy from output-\d\d\d\d/, so we exclude it to avoid duplicate.

    :param list files: A list of file in absolute path.
    :return: file list withnot the one in output-\d\d\d\d-active
    """
    files = ensure_list(files)
    active_pat = re.compile("\S*/output-(\d\d\d\d)-active/\S*")
    return [f for f in files if not active_pat.match(f)]

def filter_file(files, file_style):
    """
    Choose the file end with specified file_style

    :param list files: A list of file in absolute path.
    :param str file_style: There are few file_style you can choose:

        * .par: parameter file
        * .asc: scalar file
        * .asc: ASCII file
        * .h5: HDF5 file
        * checkpoints
    """
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



