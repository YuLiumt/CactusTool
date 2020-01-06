from .outputfile import *
import numpy as np
import json
import re


def merge_filedata(filelist):
    p = []
    for file in filelist:
        with read(file) as f:
            for dset in sorted(list(f)):
                infos = dict()
                REG = re.match('(\S+)::(\S+) it=(\d+)',dset)
                if REG:
                    infos['group'] = REG.groups()[0]
                    infos['var']   = REG.groups()[1]
                    infos['it']    = int(REG.groups()[2])
                REG = re.search('tl=(\d+)',dset); 
                if REG: 
                    infos['tl']=int(REG.groups()[0])
                REG = re.search('rl=(\d+)',dset)
                if REG: 
                    infos['rl']=int(REG.groups()[0])
                REG = re.search('c=(\d+)',dset)
                if REG: 
                    infos['c']=int(REG.groups()[0])

                subgrid = f[dset]
                try:
                    delta = subgrid.attrs['delta']
                    origin = subgrid.attrs['origin']
                    size = subgrid.shape
                    dim = len(size)
                    coord = ['x', 'y', 'z']
                    for i in range(dim) :
                        infos[coord[i]] = np.arange(0,size[(dim-1)-i])*delta[i]+origin[i]
                except:
                    print(dset)
                infos['data'] = np.array(subgrid) 
                p.append(infos)

    return p

def iteration(file):
    scalar_pat = re.compile("\S*/(output-\d\d\d\d)/\S*\.h5")
    iteration = scalar_pat.match(file)
    if iteration is not None:
        return iteration.group(1)
    else:
        return "output-0000"


class Variable:
    def __init__(self, files, var):
        self.files = files
        self.var = var
        self.iteration = {}
        for file in self.files:
            self.iteration.setdefault(iteration(file), []).append(file)

        self.Alldata = merge_filedata(self.files)
        # self.data = remove_duplicate_iters(Alldata)

    def read(self, meshgrid='HYDROBASE::press it=0 tl=0 rl=0 c=10'):
        with read(self.files[0]) as f:
            mesh = f[meshgrid]
            delta = mesh.attrs['delta']
            origin = mesh.attrs['origin']
            sizeA = mesh.shape
            tmpX = np.arange(0,sizeA[1])*delta[0]+origin[0]
            tmpY = np.arange(0,sizeA[0])*delta[1]+origin[1]

            grid = np.meshgrid(tmpX, tmpY)
            data = np.array(mesh) 
        return grid, data

    # @property
    # def grid(self):


    # @property
    # def datasets(self):

    def __str__(self):
        return "{}".format(json.dumps(self.iteration, sort_keys=True, indent=4))


def var_header(file):
    parser = re.compile(r'([^:]+)::(\S+) it=(\d+) tl=(\d+)( m=0)? rl=(\d+)( c=(\d+))?')
    vars = []
    with read(file) as f:
        for level in sorted(list(f)):
            m = parser.match(level)
            if m is not None:
                var = m.group(2)
                if var not in vars:
                    vars.append(var)
    return vars


class Griddim:
    def __init__(self, files, dim):
        self.dim = dim
        pat_fn = re.compile("\S*\.([xyz]*)\.h5(\.(gz|bz2))?$")
        self.files = [file for file in files if pat_fn.match(file).group(1) == self.dim]

        self.vars = {}
        for file in self.files:
            for var in var_header(file):
                self.vars.setdefault(var, []).append(file)

    def __getitem__(self, key):
        if key in self.vars:
            return Variable(self.vars[key], key)
        else:
            raise Exception("{} is not exist in dim {}".format(key, self.dim))
    
    def __contains__(self, key):
        return key in self.vars

    def __str__(self):
        return "Available grid data of dimension %s: \n%s\n" % (self.dim, list(self.vars.keys()))


class H5Base:
    def __init__(self, Sim):
        self.files = Sim.h5files
        self.x     = Griddim(self.files, 'x')
        self.y     = Griddim(self.files, 'y')
        self.z     = Griddim(self.files, 'z')
        self.xy    = Griddim(self.files, 'xy')
        self.xz    = Griddim(self.files, 'xz')
        self.yz    = Griddim(self.files, 'yz')
        self.xyz   = Griddim(self.files, 'xyz')

    def __str__(self):
        return "%s\n%s\n%s\n%s\n%s\n%s\n%s\n" % (self.x, self.y, self.z, self.xy, self.xz, self.yz, self.xyz)