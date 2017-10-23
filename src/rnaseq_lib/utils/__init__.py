import errno
import os

from rnaseq_lib.utils.expando import Expando

_iter_types = (list, tuple, set, frozenset)


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


def which(name, path=None):
    """
    Look for an executable file of the given name in the given list of directories,
    or the directories listed in the PATH variable of the current environment. Roughly the
    equivalent of the `which` program. Does not work on Windows.

    :type name: str
    :param name: the name of the program

    :type path: Iterable
    :param path: the directory paths to consider or None if the directories referenced in the
    PATH environment variable should be used instead

    :returns: an iterator yielding the full path to every occurrance of an executable file of the
    given name in a directory on the given path or the PATH environment variable if no path was
    passed

    >>> next( which('ls') )
    '/bin/ls'
    >>> list( which('asdalskhvxjvkjhsdasdnbmfiewwewe') )
    []
    >>> list( which('ls', path=()) )
    []
    """
    if path is None:
        path = os.environ.get('PATH')
        if path is None:
            return
        path = path.split(os.pathsep)
    for bin_dir in path:
        executable_path = os.path.join(bin_dir, name)
        if os.access(executable_path, os.X_OK):
            yield executable_path


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
        k = _key_to_attribute(k)
        if isinstance(v, dict):
            e[k] = rexpando(v)
        elif isinstance(v, _iter_types):
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
        if isinstance(v, dict):
            l.append(rexpando(v))
        elif isinstance(v, _iter_types):
            l.append(_rexpando_iter_helper(v))
        else:
            l.append(v)
    return l


def _key_to_attribute(key):
    """
    Processes key for attribute accession by replacing illegal chars with a single '_'

    :param str key: Dictionary key to process
    :return: Processed key
    :rtype: str
    """
    illegal_chars = ['!', '@', '#', '$', '%', '^', '&', '*', '(', ')', '-', '.', ',', '+', '/', '\\', ':', ';']
    for c in illegal_chars:
        key = key.replace(c, '_')
    return '_'.join(x for x in key.split('_') if x)  # Remove superfluous '_' chars
