# CarpetIOScalar
class ScalarReduction:
    """
    For gird function, We need choose which type of reduction operations. Then choose the variable we want to processes, this will pointer to another class :py:class:`Variable`.

    :param list files: A list of file in absolute path.
    :param str kind: Type of reduction operations.
    """
    def __init__(self, files, kind):
        pat_fn = re.compile("\S*\.(minimum|maximum|norm1|norm2|average)?\.asc(\.(gz|bz2))?$")
        self.kind = kind
        self.files = []
        for file in files:
            try: 
                if pat_fn.match(file).group(1) == self.kind:
                    self.files.append(file)
            except AttributeError:
                continue

    @property
    def vars(self):
        """
        All available variable in such reduction.

        :return: A dict. key is the variable name, value is corresponding files.
        """
        Vars = {}
        for file in self.files:        
            with read(file) as f:
                vars=[]
                for line in f.readlines():
                    if "# data columns: " in line:
                        vars = vars + line.split()[3:]
                        break
            assert len(vars) > 0, "{}'s header fail to identify.".format(os.path.basename(file))

            for c, name in enumerate(vars):
                Vars.setdefault(name.split(":")[1], []).append(file)

        return Vars

    def __getitem__(self, key):
        p = {}
        for k, v in self.vars.items():
            if key in k:
                p[k] = v
        if len(p) == 0:
            raise Exception("{} is not exist in reduction {}".format(key, self.kind))
        return Variable(p)

    def __contains__(self, key):
        return key in self.vars

    def __str__(self):
        if self.vars:
            return "Available %s timeseries:\n%s\n" % (str(self.kind).lower(), list(self.vars.keys()))
        else:
            return ""


class Variable:
    """
    For scalar, we don't need consider grid structure. These data may store in different files, we just simple combine them and remove duplicate data. This will done by pandas DataFrame. Sometimes you want to process a vector or a tensor. :py:class:`Variable` can also handel it.

    :param list varfiles: A dict about variable and its file, this variable may be a vector or tensor.

    * :py:attr:`Variable.Table` source data
    """
    def __init__(self, varfiles):
        self.varfiles = varfiles
        self._init_dataset()

    def _init_dataset(self):
        """
        CarpetIOScalar Dataset will store in pandas dataframe, because the ASCII data structure more like a table. This dataset will store in :py:attr:`Variable.dataset`.

        :return: DataFrame
        """
        files = []
        p = pd.DataFrame()
        for var in self.varfiles.keys():
            tem = []
            for file in self.varfiles[var]:
                filename = os.path.basename(file)
                if filename in files:
                    continue
                files.append(filename)
                data = np.loadtxt(file, comments="#")
                column = columns_asc(file)
                tem_c = pd.DataFrame(data, columns=column)
                tem.append(tem_c)
            # concat different component in same variable
            if len(tem) == 0:
                continue
            else:
                tem_p = pd.concat(tem).drop_duplicates()
            # merge different variable
            if p.empty:
                p = tem_p
            else:
                p = pd.merge(p, tem_p, how='outer', on=['time','iteration'])
 
        self.dataset = p




class CarpetIOScalar:
    """
    Thorn CarpetIOScalar provides I/O methods for outputting scalar values in ASCII format into files. This class can handle all CarpetIOScalar output. Read the dataset from files when given specified reduction operation.

    :py:class:`CarpetIOScalar` itself do nothing. You need specify the reduction operation. This will pointer to another class :py:class:`ScalarReduction`.

    :param list files: can either be a list of filenames or a single filename
    """
    def __init__(self, files):
        self.files = ensure_list(files)
        self.type = {
            'min': 'minimum',
            'max': 'maximum',
            'norm1': 'norm1',
            'norm2': 'norm2',
            'average': 'average',
            'none': None
        }

    @property
    def available_reductions(self):
        """
        available reduction operations in this dataset.

        :return: list of available reduction operations
        """
        p = []
        for type in self.type.keys():
            if str(self[type]):
                p.append(type)
        return p

    def __getattr__(self, reduction):
        """
        Specify the reduction operation by attribute

        :param str reduction: the reduction operation
        """
        assert reduction in self.type.keys(), "Does not include {} operation on scalar".format(reduction)
        return ScalarReduction(self.files, self.type[reduction])

    def __getitem__(self, reduction):
        """
        Specify the reduction operation by item

        :param str reduction: the reduction operation
        """
        assert reduction in self.type.keys(), "Does not include {} operation on scalar".format(reduction)
        return ScalarReduction(self.files, self.type[reduction])

    def __str__(self):
        return "%s%s%s%s%s%s" % (self.min, self.max, self.norm1, self.norm2, self.average, self.none)


#     @property
#     def dataset(self):
#         """
#         :return: DataFrame
#         """
#         tem = []
#         for file in self.files:
#             data = np.loadtxt(file, comments="#")
#             column = columns_asc(file)
#             tem_c = pd.DataFrame(data, columns=column)
#             tem.append(tem_c)
#         # concat different component in same variable
#         return pd.concat(tem).drop_duplicates()


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