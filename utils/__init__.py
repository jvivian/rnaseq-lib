import os

import errno


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
