import numpy as np

def l2fc(a, b):
    """
    Calculate the log2 Fold Change between two arrays, floats, or integers
    a and b cannot be, nor contain, values less than 0

    :param (int/float/np.array) a: Value or array
    :param (int/float/np.array) b: Value or array
    :return: L2FC array or value
    :rtype: (int/float/np.array)
    """
    return np.log2(a + 1) - np.log2(b + 1)
