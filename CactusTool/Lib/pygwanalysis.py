def InitModeArray(lmax):
    """Creates an array with empty modes up to lmax"""
    res = []
    for ll in range(0, lmax+1):
        res.append([])
        for mm in range(-ll, ll+1):
            res[-1].append([])
    return res