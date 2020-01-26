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
