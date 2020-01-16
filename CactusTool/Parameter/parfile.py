import os
import re


def guess_par_type(varstr, filters):
    for f in filters:
        try:
            return f(varstr)
        except ValueError:
            pass
    return str(varstr)

def par_to_str(v):
    pat = re.compile(r'^"(.*)"$', re.MULTILINE | re.DOTALL) 
    m   = re.match(pat, v)
    if m:
        return m.group(1)
    raise ValueError("Cannot convert parameter %s to string." % v)

def par_to_bool(v):
    pat_true  = re.compile(r'^"?(yes|true)"?$', re.IGNORECASE)
    pat_false = re.compile(r'^"?(no|false)"?$', re.IGNORECASE)
    if re.match(pat_true, v):
        return True
    if re.match(pat_false, v):
        return False
    raise ValueError("Cannot convert parameter to bool: %s" % v)

def par_to_varlist(v):
    s = par_to_str(v)
    pat = re.compile(r'([^\s=:"\'#!\]\[{}]+)::([^\s=:"\'#!\]\[{}]+)(\[[\t ]*[\d]+[\t ]*\])?([\t ]*\{(?:[^{}#!]|(?:\{[^{}#!]+\}))+\})?', re.MULTILINE)
    res = re.sub(pat, '', s).strip()
    if len(res) > 0:
        print(s)
        print(repr(res))
        raise ValueError("Cannot convert parameter to CACTUS variable list.")
    l = [(t.lower(),p.lower(),i,o) for t, p, i, o in pat.findall(s)]
    return set([("%s::%s%s%s" % e) for e in l])

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

def load_parfile(file):
    """
    Return:
        A dict about key value
    """
    # parameter like: 
    par_pat = re.compile(r'^[\t ]*([^\s=:"\'#!\]\[]+)::([^\s=:"\'#!\]\[]+)[\t ]*=[\t ]*([^\s=:"\'#!]+)[\t ]*(?:!|#|\n|\r\n)', re.MULTILINE)
    # parameter like: 
    vecpar_pat = re.compile(r'^[\t ]*([^\s=:"\'#!\]\[]+)::([^\s=:"\'#!\]\[]+)[\t ]*\[[\t ]*([\d]+)[\t ]*\][\t ]*=[\t ]*([^\s=:"\'#!]+)[\t ]*(?:!|#|\n|\r\n)', re.MULTILINE)
    # parameter like: Cactus::cctk_run_title = "Neutron Star"
    strpar_pat = re.compile(r'^[\t ]*([^\s=:"\'#!\]\[]+)::([^\s=:"\'#!\]\[]+)[\t ]*=[\t ]*("[^"#!]*")[\t ]*(?:!|#|\n|\r\n)', re.MULTILINE)
    # parameter like: 
    strvec_pat = re.compile(r'^[\t ]*([^\s=:"\'#!\]\[]+)::([^\s=:"\'#!\]\[]+)[\t ]*\[[\t ]*([\d]+)[\t ]*\][\t ]*=[\t ]*("[^"#!]*")[\t ]*(?:!|#|\n|\r\n)', re.MULTILINE)
    # ActiveThorn
    ath_pat = re.compile(r'^[\t ]*[A|a]ctive[T|t]horns[\t ]*=[\t ]*"([^"#]+)"[\t ]*(?:!|#|\n|\r\n)', re.MULTILINE)
    # Comment
    cmt_pat = re.compile(r'#.*')

    with open(file, 'r') as f:
        fs = f.read()
        fs = re.sub(cmt_pat, '', fs)

        athorns = ath_pat.findall(fs)
        fs      = re.sub(ath_pat, '', fs)
        strs    = strvec_pat.findall(fs)
        fs      = re.sub(strvec_pat, '', fs)
        str     = strpar_pat.findall(fs)
        fs      = re.sub(strpar_pat, '', fs)
        vector  = vecpar_pat.findall(fs)
        fs      = re.sub(vecpar_pat, '', fs)
        standard= par_pat.findall(fs)
        fs      = re.sub(par_pat, '', fs).strip()

        if len(fs)>0:
            prin("Warning: unparsed parfile content")
            print(fs)

    known_varlists = set([
        ('iobasic','outinfo_vars'),
        ('ioscalar','outscalar_vars'),
        ('ioascii','out0d_vars'),
        ('ioascii','out1d_vars'),
        ('iohdf5','out1d_vars'),
        ('iohdf5','out2d_vars'),
        ('iohdf5','out3d_vars'),
        ('iohdf5','out_vars'),
        ('carpetiobasic','outinfo_vars'),
        ('carpetioscalar','outscalar_vars'),
        ('carpetioascii','out0d_vars'),
        ('carpetioascii','out1d_vars'),
        ('carpetiohdf5','out1d_vars'),
        ('carpetiohdf5','out2d_vars'),
        ('carpetiohdf5','out3d_vars'),
        ('carpetiohdf5','out_vars'),
        ('dissipation', 'vars'),
        ('nanchecker', 'check_vars'),
        ('summationbyparts', 'vars')
    ])
 
    filters = [par_to_bool, int, float, par_to_str]
    parfilt = lambda s: guess_par_type(s, filters)

    p = dict()

    for thorn, param, value in str:
        if (thorn.lower(), param.lower()) in known_varlists:
            add_two_key(p, thorn, param, par_to_varlist(value))
        else:
            add_two_key(p, thorn, param, parfilt(value))

    for thorn, param, value in standard:
        add_two_key(p, thorn, param, parfilt(value))
        
    vector_pattern = vector + strs
    for thorn, param, idx, value in vector_pattern:
        add_three_key(p, thorn, param, idx, parfilt(value))

    p['ActiveThorns'] = []
    for thlist in athorns:
        for t in thlist.split():
            p['ActiveThorns'].append(t)
            if t not in p.keys():
                p.update({t: None})

    return p


class Thorn:
    """
    This class represents the parameters of a given Cactus thorn.
    """
    def __init__(self, parameter):
        self.parameter = parameter

    def __str__(self):
        return "%s" % (self.parameter)


class ParFile:
    def __init__(self, parfiles):
        # Use the parfile which is first created
        par_create_time = [os.path.getctime(par) for par in parfiles]
        self.parfile = parfiles[par_create_time.index(min(par_create_time))]
    
        self.dict = load_parfile(self.parfile)

    @property
    def activethorn(self):
        return self.dict["ActiveThorns"]

    def __getattr__(self, name):
        if name in self.dict["ActiveThorns"]:
            return Thorn(self.dict[name])
        else:
            raise ValueError("Thorn %s not in ActiveThorn. Please check the letter case." % name)

    def __str__(self):
        return "%s:\n%s" % (os.path.basename(self.parfile), sorted(self.dict.keys()))