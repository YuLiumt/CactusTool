import numpy as np

def is_empty(obj):
    """
    Make sure `obj` isn't empty.
    """
    if not obj:
        print("Please check why {} is None".format(obj))

def ensure_list(obj):
    """
    This function ensures that `obj` is a list. Typically used to convert a string to a list.
    """
    if obj is None:
        return [obj]
    if not isinstance(obj, list):
        return [obj]
    return obj

def ensure_numpy_array(obj):
    """
    This function ensures that `obj` is a numpy array. Typically used to convert scalar, list or tuple argument passed to functions using Cython.
    """
    if isinstance(obj, np.ndarray):
        if obj.shape == ():
            return np.array([obj])
        # We cast to ndarray to catch ndarray subclasses
        return np.array(obj)
    elif isinstance(obj, (list, tuple)):
        return np.asarray(obj)
    else:
        return np.asarray([obj])