import errno
import os
from collections import OrderedDict
from collections import defaultdict

from rnaseq_lib.utils.expando import Expando

dict_types = (dict, OrderedDict, defaultdict)
iter_types = (list, tuple, set, frozenset)


def mkdir_p(path):
    """
    Creates directory unless it already exists

    :param str path: Path of directory to make
    """
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def flatten(x):
    """
    Flattens a nested array into a single list

    :param list x: The nested list/tuple to be flattened
    """
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result


def merge_dicts(x, y):
    """
    Given two dicts, merge them into a new dict as a shallow copy

    param dict x: first dictionary
    param dict y: second dictionary
    """
    z = x.copy()
    z.update(y)
    return z


def rexpando(d):
    """
    Recursive Expando!

    Recursively iterate through a nested dict / list object
    to convert all dictionaries to Expando objects

    :param dict d: Dictionary to convert to nested Expando objects
    :return: Converted dictionary
    :rtype: Expando
    """
    e = Expando()
    for k, v in d.iteritems():
        if any(isinstance(v, t) for t in dict_types):
            e[k] = rexpando(v)
        elif any(isinstance(v, t) for t in iter_types):
            e[k] = _rexpando_iter_helper(v)
        else:
            e[k] = v
    return e


def _rexpando_iter_helper(input_iter):
    """
    Recursively handle iterables for rexpando

    :param iter input_iter: Iterable to process
    :return: Processed iterable
    :rtype: list
    """
    l = []
    for v in input_iter:
        if any(isinstance(v, t) for t in dict_types):
            l.append(rexpando(v))
        elif any(isinstance(v, t) for t in iter_types):
            l.append(_rexpando_iter_helper(v))
        else:
            l.append(v)
    return l
