from .check import *
import pandas as pd
import h5py
import copy
import re
import os


################### Folder ###################
def fetch_all_file(path):
    """
    Fetch all files under the path.

    :param str path: directory in absolute path
    :return: file list
    """
    assert os.path.exists(path), "{} doesn't exist.".format(path)

    filelist = []
    for root, dirs, files in os.walk(path):
        for file in files:
            filelist.append(os.path.join(root, file))

    return filelist

def rm_dir_in_filelist(files, dir):
    """
    Deletes files under the specified directory.

    :param list files: file list in absolute path
    :param str dir: The directory you want to remove
    :return: file list
    """
    files = ensure_list(files)
    dir_pat = re.compile("\S*/" + dir + "/\S*")
    return [f for f in files if not dir_pat.match(f)]

def filter_file(files, file_style):
    """
    Choose the file end with specified file_style

    :param list files: file list in absolute path.
    :param str file_style: There are few file_style you can choose:

        * parfile: parameter file
        * scalar: Scalar file
        * ascii: ASCII file
        * hdf5: HDF5 file
        * debug: NaNmask
    
    :return: file list
    """
    files = ensure_list(files)
    re_pat = {
        "parfile": "\S*\.par",
        "scalar": "\S*\.(minimum|maximum|norm1|norm2|norm_inf|average)?\.asc(\.(gz|bz2))?$",
        "ascii": "\S*\.[xyz]+\.asc(\.(gz|bz2))?$",
        "hdf5": "\S*\.[xyz]+(\.file_\d+)?\.h5(\.(gz|bz2))?$",
        "array": "\S*\.asc$",
        "debug": "\S*NaNmask\.\S*\.h5"
    }
    return [f for f in files if re.compile(re_pat[file_style]).match(f)]


def read(file):
    """
    Carpet output file type is completely different. This function provides different way to open it.

    :param str file: file in absolute path

    >>> with read(file) as f:
    """
    assert os.path.exists(file), "{} doesn't exist in your local computer".format(file)
    if file.endswith('.h5'):
        f = h5py.File(file, 'r')
    elif file.endswith(('.par', '.asc', '.txt')):
        f = open(file)
    else:
        raise RuntimeError("CactusTool can't handle this *{} type of file".format(os.path.splitext(file)[1]))

    return f



# def dataset(file, slice, dim):
#     p = dict()

#     # with read(file) as f:
#     #     if slice not in sorted(list(f)):
#     #         raise RuntimeError("%s not in %s" % (slice, file))
#     #     dset = f[slice]
#     #     for item in list(dset.attrs):
#     #         p[item] = dset.attrs.get(item, None)
#     #     p['data'] = np.array(dset)  

#     with read(file) as f:
#         if slice not in sorted(list(f)):
#             raise RuntimeError("%s not in %s" % (slice, file))
#         mesh = f[slice]
#         p['level'] = mesh.attrs.get('level', None)
#         origin = mesh.attrs.get('origin', None)
#         delta = mesh.attrs.get('delta', None)
#         p['origin'] = origin - delta/2
#         data = np.array(mesh)
#         size = mesh.shape
#         n = len(dim)
#         for i in range(n):
#             p[dim[i]] = np.arange(0, size[(n-1)-i]+1)*delta[i] + origin[i]
#         # p['right_edge'] = [size[(dim-1)-i]*delta[i]+origin[i] for i in range(dim)]
#         # p['dimensions'] = size
#         p['data'] = data

#     return p

# def dataframe_h5(file, axis):
#     """
#     HDF5 to pandas DataFrame
#     """
#     assert axis in ['x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz'], "Does not include {} axis on {}".format(dim, file)
#     pattern = re.compile(r'([^:]+)::(\S+) it=(\d+) tl=(\d+)( m=(\d+))? rl=(\d+)( c=(\d+))?')
#     with read(file) as f:
#         tem = []
#         for item in list(f):
#             header = pattern.match(item)
#             if header is None:
#                 continue
#             dset = f[item]
#             var = header.group(2)
#             data = np.array(dset)
#             origin = dset.attrs.get('origin', None)
#             delta = dset.attrs.get('delta', None)
#             size = data.shape
#             dim = len(axis)
#             coord = tuple(np.arange(0,size[(dim-1)-i])*delta[i]+origin[i] for i in range(dim))
#             grid = np.meshgrid(*coord)
#             tem_c = pd.DataFrame()
#             for i in range(dim):
#                 tem_c[axis[i]] = grid[i].flatten()
#             tem_c[var] = data.flatten()
#             tem_c['rl'] = int(header.group(7))
#             tem_c['it'] = int(header.group(3))
#             tem_c['time'] = dset.attrs.get('time', None)
#             if header.group(9) != None:
#                 tem_c['c'] = int(header.group(9))    
#             tem.append(tem_c)

#     return pd.concat(tem)

##### HDF5 #####
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
                # 'file': file,
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

def select_header_h5(header, var=-1, it=-1, ml=-1, rl=-1, c=-1):
    """
    Select specified header in return of :py:func:`header_h5`.

    :param str var: -1: all variable.
    :param int it: -1: all iteration number.
    :param int ml: -1: all map.
    :param int rl: -1: all refinement level.
    :param int c: -1: all component.
    :return: a list of header
    """
    p = copy.deepcopy(header)
    for f, hs in header.items():
        for h, item in hs.items(): 
            if var != -1 and item['varname'] != var:
                del p[f][h]
                continue
            if it != -1 and item['iteration'] != it:
                del p[f][h]
                continue
            if ml != -1 and item['ml'] != ml:
                del p[f][h]
                continue
            if rl != -1 and item['rl'] != rl:
                del p[f][h]
                continue
            if c != -1 and 'c' in item and item['c'] != c:
                del p[f][h]
                continue
        if not p[f]:
            del p[f]

    return p

def dataset_h5(header):
    """
    Return a dictionary about dataset with in slice.

    :param str slice: header inf
    :return: a dictionary about dataset
    """
    p = dict()
    
    for file, hs in header.items():
        with read(file) as f: 
            for h, v in hs.items():
                mesh = f[h]
                time = round(mesh.attrs.get('time', None), 6)
                if time not in p:
                    p[time] = {}
                ml = v['ml']
                if ml not in p[time]:
                    p[time][ml] = {}
                rlevel = v['rl']
                if rlevel not in p[time][ml]:
                    p[time][ml][rlevel] = {}
                if 'c' in v:
                    component = v['c']
                else:
                    component = 0
                if component not in p[time][ml][rlevel]:
                    p[time][ml][rlevel][component] = {}
                p[time][ml][rlevel][component].update({'origin': mesh.attrs.get('origin', None)})
                p[time][ml][rlevel][component].update({'delta': mesh.attrs.get('delta', None)})
                p[time][ml][rlevel][component].update({'ghostzones': mesh.attrs.get('ghostzones', None)})
                p[time][ml][rlevel][component].update({'data': np.array(mesh)})
    
    return p

##### ASCII #####
def columns_asc(file):
    """
    Extract ASCII file header information.

    :param str file: file in absolute path
    :return: The columns in given file.
    """
    with read(file) as f:
        for line in f.readlines(2000):
            if "# 1:iteration 2:time 3:data" in line:
                columns = line.split()[1:]
            if "# column format:" in line:
                columns = line.split()[3:]
            if "# data columns:" in line:
                del columns[-1]
                columns += line.split()[3:]
                break
    assert len(columns) > 0, "File: {} Header fail to identify.".format(file)
    return [name.split(":")[1] for c, name in enumerate(columns)]


def dataset_asc(files, usecols):
    p = np.empty((len(usecols), 0))
    for f in files:
        tem = np.loadtxt(f, comments="#", usecols=usecols, unpack=True)
        p = np.append(p, tem, axis=1)
    return np.unique(p, axis=1)