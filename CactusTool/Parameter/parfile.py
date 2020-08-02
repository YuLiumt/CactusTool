from loguru import logger
import os
import re


def cactus_format(value):
    if isinstance(value, list):
        return '"\n    '+'\n    '.join(value)+ '\n"'
    elif isinstance(value, bool):
        return '"yes"' if value else '"no"'
    elif isinstance(value, str):
        return '"' + value + '"'
    else:
        
        return value

def par_to_bool(v):
    pat_true  = re.compile(r'^"?(yes|true)"?$', re.IGNORECASE) 
    pat_false = re.compile(r'^"?(no|false)"?$', re.IGNORECASE) 
    if re.match(pat_true, v): return True
    if re.match(pat_false, v): return False
    raise ValueError("Cannot convert parameter %s to bool." % v)

def par_to_str(v):
    pat = re.compile(r'^"(.*)"$', re.MULTILINE | re.DOTALL) 
    m   = re.match(pat, v)
    if m:
        return m.group(1)
    raise ValueError("Cannot convert parameter %s to string." % v)

def par_format(value):
    for f in [par_to_bool, int, float, par_to_str]:
        try:
            return f(value)
        except ValueError:
            pass
    return value


class Thorn:
    """
    This class represents the parameters of a given Cactus thorn.
    """
    def __init__(self, thorn):
        self.params = {}
        self.thorn = thorn.lower()

    def __contains__(self, param):
        return str(param).lower() in self.params

    def __iter__(self):
        return iter(self.params.items())

    def __getitem__(self, key):
        key = key.lower()
        if key in self:
            return self.params[key]  
        else:
            raise KeyError("Parameter %s not set." % key)

    def __setitem__(self, key, value):
        self.params[key.lower()] = value

    # def __getattr__(self, key):
    #     return self[key]

    def __str__(self):
        p = []
        for par, val in self.params.items():
            key = "%s::%s" % (self.thorn, par.lower())
            if val == '$parfile':
                parline = "%s = %s" % (key, val)
            else:
                parline = "%s = %s" % (key, cactus_format(val))
            p.append(parline)
        return "\n".join(p)

class ParFile:

    def __init__(self, parfile):
        """
        Load parfile.

        :param str file: file in absolute path
        """
        assert os.path.exists(parfile), "{} doesn't exist.".format(parfile)

        self.parfile = parfile
        self.activethorns = set()
        self.thorns = {}
        self._load_parfile()
    
    def _load_parfile(self):
        """
        Initialize.
        """
        # Comment
        cmt_pat = re.compile(r'#.*')
        # ActiveThorn
        ath_pat = re.compile(r'^[\t ]*[A|a]ctive[T|t]horns[\t ]*=[\t ]*"([^"#]+)"[\t ]*(?:!|#|\n|\r\n)', re.MULTILINE)
        # parameter like: AHFinderDirect::track_origin_source_x[1] = "PunctureTracker::pt_loc_x[0]"
        strvec_pat = re.compile(r'^[\t ]*([^\s=:"\'#!\]\[]+)::([^\s=:"\'#!\]\[]+)[\t ]*\[[\t ]*([\d]+)[\t ]*\][\t ]*=[\t ]*("[^"#!]*")[\t ]*(?:!|#|\n|\r\n)', re.MULTILINE)
        # parameter like: CarpetRegrid2::radius_1[1] = 240.0
        vecpar_pat = re.compile(r'^[\t ]*([^\s=:"\'#!\]\[]+)::([^\s=:"\'#!\]\[]+)[\t ]*\[[\t ]*([\d]+)[\t ]*\][\t ]*=[\t ]*([^\s=:"\'#!]+)[\t ]*(?:!|#|\n|\r\n)', re.MULTILINE)
        # parameter like: Cactus::cctk_run_title = "Neutron Star"
        strpar_pat = re.compile(r'^[\t ]*([^\s=:"\'#!\]\[]+)::([^\s=:"\'#!\]\[]+)[\t ]*=[\t ]*("[^"#!]*")[\t ]*(?:!|#|\n|\r\n)', re.MULTILINE)
        # parameter like: cactus::cctk_final_time = 2500.0
        par_pat = re.compile(r'^[\t ]*([^\s=:"\'#!\]\[]+)::([^\s=:"\'#!\]\[]+)[\t ]*=[\t ]*([^\s=:"\'#!]+)[\t ]*(?:!|#|\n|\r\n)', re.MULTILINE)

        with open(self.parfile, 'r') as f:
            fs      = f.read()
            fs      = re.sub(cmt_pat, '', fs)
            athorns = ath_pat.findall(fs)
            fs      = re.sub(ath_pat, '', fs)
            pstr    = strpar_pat.findall(fs)
            fs      = re.sub(strpar_pat, '', fs)
            strvec  = strvec_pat.findall(fs)
            fs      = re.sub(strvec_pat, '', fs)
            numvec  = vecpar_pat.findall(fs)
            fs      = re.sub(vecpar_pat, '', fs)
            standard= par_pat.findall(fs)
            fs      = re.sub(par_pat, '', fs).strip()
            if len(fs)>0:
                logger.warning("unparsed parfile content:\n{}", fs)

        for thlist in athorns:
            self.activate_thorns(thlist.split())

        known_varlists = set([
            ('iobasic','outinfo_vars'),
            ('ioscalar','outscalar_vars'),
            ('ioascii','out0d_vars'),
            ('ioascii','out1d_vars'),
            ('ioascii','out2d_vars'),
            ('ioascii','out3d_vars'),
            ('iohdf5','out1d_vars'),
            ('iohdf5','out2d_vars'),
            ('iohdf5','out3d_vars'),
            ('iohdf5','out_vars'),
            ('carpetiobasic','outinfo_vars'),
            ('carpetioscalar','outscalar_vars'),
            ('carpetioascii','out0d_vars'),
            ('carpetioascii','out1d_vars'),
            ('carpetioascii','out2d_vars'),
            ('carpetioascii','out3d_vars'),
            ('carpetiohdf5','out1d_vars'),
            ('carpetiohdf5','out2d_vars'),
            ('carpetiohdf5','out3d_vars'),
            ('carpetiohdf5','out_vars'),
            ('dissipation', 'vars'),
            ('nanchecker', 'check_vars'),
            ('summationbyparts', 'vars'),
            ('hmns', 'out1d_array'),
            ('hmns', 'out2d_array'),
        ])

        for thorn, param, value in pstr:
            if (thorn.lower(), param.lower()) in known_varlists:
                pat = re.compile(r'([^\s=:"\'#!\]\[{}]+)::([^\s=:"\'#!\]\[{}]+)(\[[\t ]*[\d]+[\t ]*\])?([\t ]*\{(?:[^{}#!]|(?:\{[^{}#!]+\}))+\})?', re.MULTILINE)
                l = [(t.lower(),p.lower(),i,o) for t,p,i,o in pat.findall(value)]
                self[thorn][param] = [("%s::%s%s%s" % e) for e in l]
            else:
                self[thorn][param] = par_format(value)

        for thorn, param, idx, value in strvec:
            self[thorn]["{}[{}]".format(param, idx)] = par_format(value)

        for thorn, param, value in standard:
            self[thorn][param] = par_format(value)

        for thorn, param, idx, value in numvec:
            self[thorn]["{}[{}]".format(param, idx)] = par_format(value)

    def activate_thorns(self, thorn):
        if isinstance(thorn, list):
            th = [str(t).lower() for t in thorn]
            self.activethorns.update(th)
        else:
            th = str(thorn).lower()
            self.activethorns.add(th)

    def __getitem__(self, thorn):
        th = str(thorn).lower()
        if th in self.thorns:
            return self.thorns[th]
        else:
            self.thorns[th] = Thorn(th)
            return self.thorns[th]

    def __contains__(self, thorn):
        return str(thorn).lower() in self.thorns

    def __iter__(self):
        return iter(self.thorns)

    def __getattr__(self, thorn):
        return self[thorn]

    def __str__(self):
        output = 'ActiveThorns = "\n    ' + "\n    ".join(sorted(list(self.activethorns))) +'\n"\n\n'
        for thorn in self:
            output += "\n########## {} ##########\n".format(thorn)
            output += str(self[thorn]) + '\n'
        return output

    def save(self, filename):
        assert filename.endswith('.par'), "file name '{}' should end with '.par'.".format(filename)
        with open(filename, 'w') as f:
            f.write(str(self))
