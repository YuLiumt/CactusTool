"""
Cactus dataset main produced by `Carpet <https://carpetcode.org>`_, which is an adaptive mesh refinement and multi-patch driver for the Cactus. This modular processes output of CarpetIOHDF5 or CarpetIOASCII.
"""
from ..funcs.file import header_h5, select_header_h5, read, columns_asc
from ..funcs.json import Format, add_two_key
from ..funcs.log import logger
from ..Visualize.plot2D import pcolormesh
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import copy
import re
import os

class CarpetGF:
    """
    A class representing grid function.
    """
    def __init__(self, files, dim, ftype):
        self.files = files
        self.dim = dim
        self.ftype = ftype
        self.fname = set()
        pat = re.compile('^([a-zA-Z0-9_-]+)(\.\d+)?\.([xyz]*)(\.file_\d+)?\.(asc|h5)$')
        for f in files:
            mp = pat.match(os.path.basename(f))
            if mp is not None:
                thorn_var = mp.group(1)
                self.fname.add(thorn_var)
    
    def __getitem__(self, key):
        """
        Read the dataset when given specified file name.

        :param str key: file name.
        """
        assert key in self.fname, "{} is not exist".format(key)
        fileList = [f for f in self.files if key+'.' in f]
        if self.ftype == 'h5':
            return CarpetIOHDF5(fileList, self.dim)
        elif self.ftype == 'asc':
            assert len(self.dim) == 1, "We do not recommend using ASCII to draw high-dimensional graphs"
            return CarpetIOASCII(fileList, self.dim)

    def __contains__(self, key):
        return key in self.fname



class CarpetIOHDF5:
    """
    Thorn CarpetIOHDF5 provides I/O methods for outputting gird function in HDF5 format into files. This class can handle all CarpetIOHDF5 output. Read the dataset from files and return all useful information stored in these files. We recommend use CarpetIOHDF5 to output all 2-D and 3-D grid function.
    """
    
    def __init__(self, files, dim):
        """
        :param list files: can either be a list of filenames or a single 
        :param str dim: dim
        """
        self.dim = dim
        self.header = {}
        for file in files:
            self.header[file] = header_h5(file)

    @property
    def vars(self):
        vars = set()
        for file in self.header:
            for item in self.header[file]:
                vars.add(self.header[file][item]['varname'])
        return list(vars)

    @property
    def it(self):
        iteration = set()
        for file in self.header:
            for item in self.header[file]:
                iteration.add(int(self.header[file][item]['iteration']))
        return sorted(list(iteration))

    def dsets(self, var, it=0, rl=-1, c=-1):
        p = {}
        for file in self.header:
            with read(file) as f: 
                for header in select_header_h5(self.header[file], var, it=it, rl=rl, c=c):
                    rl = self.header[file][header]['rl']
                    c = self.header[file][header]['c']
                    dset = {}
                    mesh = f[header]
                    dset['origin'] = mesh.attrs.get('origin', None)
                    dset['delta'] = mesh.attrs.get('delta', None)
                    dset['ghostzones'] = mesh.attrs.get('cctk_nghostzones', None)
                    dset['time'] = mesh.attrs.get('time', None)
                    dset['data'] = np.array(mesh) 
                    add_two_key(p, rl, c, dset)

        return p

    def animation(self, var, axlim=None, reflevel=-1, unit='cactus', saveto=None, cmap="viridis", interval=100, **kwargs):
        """
        plot 2D Carpet data, then animate it.

        :param int interval: Pause between frames in ms
        """
        pass
        # fig, ax = plt.subplots()

        # if axlim is not None:
        #     ax.set_xlim([-axlim, axlim])
        #     ax.set_ylim([-axlim, axlim])

        # it = self.it
        # ims = pcolormesh(ax, self.dsets(var, it=it[0]))

        # def animate(n):
        #     now = it[n]
        #     ims = pcolormesh(ax, self.dsets(var, it=now))
        #     return ims
        #     line.set_data([], [])
        #     dsets = self.dsets(var, it=now)
        #     for rl in sorted(dsets):
        #         if reflevel != -1 and rl != reflevel:
        #             continue
        #         for c in sorted(dsets[rl]):
        #             coord = dsets[rl][c]['coord']
        #             data = dsets[rl][c]['data']
        #             time = dsets[rl][c]['time']
        #             # Add the data to the line object
        #             line.set_xdata(np.append(line.get_xdata(), coord))
        #             line.set_ydata(np.append(line.get_ydata(), data))
        #     # Sort points
        #     indices = np.argsort(line.get_xdata())
        #     line.set_xdata(line.get_xdata()[indices])
        #     line.set_ydata(line.get_ydata()[indices])
        #     # Adjust the axes
        #     ax.relim()
        #     ax.autoscale_view()
        #     if unit == 'cactus':
        #         ax.set_xlabel('Coord [M]')
        #     ax.set_title("Time: {}".format(time))
        #     return line,

        # anim = animation.FuncAnimation(fig, animate, frames=len(it), interval=interval, blit=True, repeat=False)
        
        # if saveto is not None:
        #     # writer = animation.FFMpegWriter(fps=15, metadata=dict(artist='Me'), bitrate=1800)
        #     anim.save(saveto)
        # return anim

class CarpetIOASCII:
    """
    Thorn CarpetIOASCII provides I/O methods for outputting gird function in ASCII format into files. This class can handle all CarpetIOASCII output. Read the dataset from files and return all useful information stored in these files.

    :param str files: a single filename
    """
    def __init__(self, file, dim):
        self.dim = dim
        self.columns = columns_asc(file)
        self.vars = self.columns[12:]
        self._dsets = np.loadtxt(file, comments="#")

    @property
    def it(self):
        it = self._dsets.T[0].astype(np.int)
        return sorted(list(set(it)))

    def dsets(self, var='data', it=0):
        p = {}
        dsets = self._dsets[self._dsets[:, 0] == it]
        for rl in set(dsets[:, 2]):
            dset = dsets[dsets[:, 2] == rl]
            for c in set(dset[:, 3]):
                data = dset[dset[:, 3] == c].T
                p_tem = {}
                if self.vars[-1] == 'data':
                    p_tem['data'] = data[12]
                else:
                    assert var in self.vars, "Pick one from {}, because one file per group.".format(self.vars)
                    p_tem['data'] = data[self.columns.index(var)]
                p_tem['coord'] = data[self.columns.index(self.dim)]
                p_tem['time'] = data[8][0]
                add_two_key(p, rl, c, p_tem)

        return p

    def animation(self, var='data', reflevel=-1, unit='cactus', saveto=None, **kwargs):
        """
        plot 1D Carpet data, then animate it.
        """
        fig, ax = plt.subplots()
        it = self.it
        line = ax.plot([], [], **kwargs)[0]

        def animate(n):
            now = it[n]
            line.set_data([], [])
            dsets = self.dsets(var, it=now)
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
            return line,

        anim = animation.FuncAnimation(fig, animate, frames=len(it), interval=1000, blit=True, repeat=False)
        
        if saveto is not None:
            # writer = animation.FFMpegWriter(fps=15, metadata=dict(artist='Me'), bitrate=1800)
            anim.save(saveto)
        return anim



# origin = mesh.attrs.get('origin', None)
# delta = mesh.attrs.get('delta', None)
# data = np.array(mesh)
# size = mesh.shape
# n = len(self.dim)
# for i in range(n):
#     dset[self.dim[i]] = np.arange(0, size[(n-1)-i]+1)*delta[i] + origin[i] - delta[i]/2
# dset['data'] = data        


#     @property
#     def rl(self):
#         rl = set()
#         for file in self.header:
#             for item in self.header[file]:
#                 rl.add(int(self.header[file][item]['rl']))
#         return sorted(list(rl))

#     def temporary(self, var, it=0, rl=-1):
#         dsets = []
#         if rl == -1:
#             rl = self.rl
#         for level in list(rl):
#             for file in self.header:
#                 with read(file) as f: 
#                     headers = select_header_h5(self.header[file], var, it=it, rl=level)
            
#                     for slice in headers:
#                         p = dict()
#                         mesh = f[slice]
#                         p['level'] = mesh.attrs.get('level', None)
#                         origin = mesh.attrs.get('origin', None)
#                         delta = mesh.attrs.get('delta', None)
#                         p['origin'] = origin
#                         data = np.array(mesh)
#                         size = mesh.shape
#                         n = len(self.dim)
#                         for i in range(n):
#                             p[self.dim[i]] = np.arange(0, size[(n-1)-i]+1)*delta[i] + origin[i] - delta[i]/2
#                         p['data'] = data
#                         dsets.append(p)
#         return AMRGrid(dsets, self.dim, var)

    # def eval(self, expr, it):
    #     var = None
    #     if ('gxx' or 'gxy' or 'gxz' or'gyy' or'gyz' or'gzz') in expr:
    #         admbase_metric = self['admbase-metric']
    #         gxx = admbase_metric.temporary('gxx', it=it)
    #         gxy = admbase_metric.temporary('gxy', it=it)
    #         gxz = admbase_metric.temporary('gxz', it=it)
    #         gyy = admbase_metric.temporary('gyy', it=it)
    #         gyz = admbase_metric.temporary('gyz', it=it)
    #         gzz = admbase_metric.temporary('gzz', it=it)
    #         expr = expr.replace('gxx', "gxx[rl][c]['data']").replace('gxy', "gxy[rl][c]['data']").replace('gxz', "gxz[rl][c]['data']").replace('gyy', "gyy[rl][c]['data']").replace('gyz', "gyz[rl][c]['data']").replace('gzz', "gzz[rl][c]['data']")
    #         var = 'gxx'
    #     if 'alp' in expr:
    #         admbase_lapse = self['admbase-lapse']
    #         alp = admbase_lapse.temporary('alp', it=it)
    #         expr = expr.replace('alp', "alp[rl][c]['data']")
    #         var = 'alp'
    #     if ('betax' or 'betay' or 'betax') in expr:
    #         admbase_shift = self['admbase-shift']
    #         betax = admbase_shift.temporary('betax', it=it)
    #         betay = admbase_shift.temporary('betay', it=it)
    #         betaz = admbase_shift.temporary('betaz', it=it)
    #         expr = expr.replace('betax', "betax[rl][c]['data']").replace('betay', "betay[rl][c]['data']").replace('betaz', "betaz[rl][c]['data']")
    #         var = 'betax'
    #     if 'rho' in expr:
    #         hydrobase_rho = self['hydrobase-rho']
    #         rho = hydrobase_rho.temporary('rho', it=it)
    #         expr = expr.replace('rho', "rho[rl][c]['data']")
    #         var = 'rho'
    #     if 'press' in expr:
    #         hydrobase_press = self['hydrobase-press']
    #         press = hydrobase_press.temporary('press', it=it)
    #         expr = expr.replace('press', "press[rl][c]['data']")
    #         var = 'press'
    #     if 'eps' in expr:
    #         hydrobase_eps = self['hydrobase-eps']
    #         eps = hydrobase_eps.temporary('eps', it=it)
    #         expr = expr.replace('eps', "eps[rl][c]['data']")
    #         var = 'eps'
    #     if 'w_lorentz' in expr:
    #         hydrobase_w_lorentz = self['hydrobase-w_lorentz']
    #         w_lorentz = hydrobase_w_lorentz.temporary('w_lorentz', it=it)
    #         expr = expr.replace('w_lorentz', "w_lorentz[rl][c]['data']")
    #         var = 'w_lorentz'
    #     if ('velx' or 'vely' or 'velz') in expr:
    #         hydrobase_vel = self['hydrobase-vel']
    #         velx = hydrobase_vel.temporary('vel[0]', it=it)
    #         vely = hydrobase_vel.temporary('vel[1]', it=it)
    #         velz = hydrobase_vel.temporary('vel[2]', it=it)
    #         expr = expr.replace('velx', "velx[rl][c]['data']").replace('vely', "vely[rl][c]['data']").replace('velz', "velz[rl][c]['data']")
    #         var = 'velx'

    #     p = copy.deepcopy(locals()[var])
    #     for rl in p:
    #         for c in p[rl]:
    #             p[rl][c]['data'] = eval(expr)
    #     return p