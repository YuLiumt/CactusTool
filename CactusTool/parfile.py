import re

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
    
def par_to_str(v):
    pat = re.compile(r'^"(.*)"$', re.MULTILINE | re.DOTALL) 
    m   = re.match(pat, v)
    if m:
        return m.group(1)
    raise ValueError("Cannot convert parameter %s to string." % v)

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

def guess_par_type(varstr, filters, default):
    for f in filters:
        try:
            return f(varstr)
        except ValueError:
            pass
    return default(varstr)

def par_to_bool(v):
    pat_true  = re.compile(r'^"?(yes|true)"?$', re.IGNORECASE)
    pat_false = re.compile(r'^"?(no|false)"?$', re.IGNORECASE)
    if re.match(pat_true, v):
        return True
    if re.match(pat_false, v):
        return False
    raise ValueError("Cannot convert parameter to bool: %s" % v)

def load_parfile(file):
    par_pat = re.compile(r'^[\t ]*([^\s=:"\'#!\]\[]+)::([^\s=:"\'#!\]\[]+)[\t ]*=[\t ]*([^\s=:"\'#!]+)[\t ]*(?:!|#|\n|\r\n)', re.MULTILINE)
    vecpar_pat = re.compile(r'^[\t ]*([^\s=:"\'#!\]\[]+)::([^\s=:"\'#!\]\[]+)[\t ]*\[[\t ]*([\d]+)[\t ]*\][\t ]*=[\t ]*([^\s=:"\'#!]+)[\t ]*(?:!|#|\n|\r\n)', re.MULTILINE)
    strpar_pat = re.compile(r'^[\t ]*([^\s=:"\'#!\]\[]+)::([^\s=:"\'#!\]\[]+)[\t ]*=[\t ]*("[^"#!]*")[\t ]*(?:!|#|\n|\r\n)', re.MULTILINE)
    strvec_pat = re.compile(r'^[\t ]*([^\s=:"\'#!\]\[]+)::([^\s=:"\'#!\]\[]+)[\t ]*\[[\t ]*([\d]+)[\t ]*\][\t ]*=[\t ]*("[^"#!]*")[\t ]*(?:!|#|\n|\r\n)', re.MULTILINE)
    ath_pat = re.compile(r'^[\t ]*[A|a]ctive[T|t]horns[\t ]*=[\t ]*"([^"#]+)"[\t ]*(?:!|#|\n|\r\n)', re.MULTILINE)
    cmt_pat = re.compile(r'#.*')

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

    with open(file, 'r') as f:
        fs = f.read()
        fs = re.sub(cmt_pat, '', fs)
        athorns = ath_pat.findall(fs)
        fs      = re.sub(ath_pat, '', fs)
        pstrvec = strvec_pat.findall(fs)
        fs      = re.sub(strvec_pat, '', fs)
        pstr    = strpar_pat.findall(fs)
        fs      = re.sub(strpar_pat, '', fs)
        pvec    = vecpar_pat.findall(fs)
        fs      = re.sub(vecpar_pat, '', fs)
        pstd    = par_pat.findall(fs)
        fs      = re.sub(par_pat, '', fs).strip()
        if len(fs)>0:
            prin("Warning: unparsed parfile content")
            print(fs)
    
    filters = [par_to_bool, int, float, par_to_str]
    parfilt = lambda s: guess_par_type(s, filters, str)

    p = dict()
    for thorn, param, value in pstr:
        if (thorn.lower(), param.lower()) in known_varlists:
            add_two_key(p, thorn, param, par_to_varlist(value))
            # print(par_to_varlist(value))
        else:
            add_two_key(p, thorn, param, parfilt(value))

    for thorn, param, value in pstd:
        add_two_key(p, thorn, param, parfilt(value))
    
    vpa = pvec + pstrvec
    for thorn, param, idx, value in vpa:
        add_three_key(p, thorn, param, idx, parfilt(value))

    p['ActiveThorns'] = []
    for thlist in athorns:
        for t in thlist.split():
            p['ActiveThorns'].append(t)
            if t not in p.keys():
                p.update({t: None})

    return p

