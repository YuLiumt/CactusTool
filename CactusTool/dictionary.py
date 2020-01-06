import json


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