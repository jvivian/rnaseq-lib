import numpy as np


# Outlier
def iqr_bounds(ys):
    """
    Return upper and lower bound for an array of values

    Lower bound: Q1 - (IQR * 1.5)
    Upper bound: Q3 + (IQR * 1.5)

    :param list ys: List of values to calculate IQR
    :return: Upper and lower bound
    :rtype: tuple(float, float)
    """
    quartile_1, quartile_3 = np.percentile(ys, [25, 75])
    iqr = quartile_3 - quartile_1
    lower_bound = quartile_1 - (iqr * 1.5)
    upper_bound = quartile_3 + (iqr * 1.5)
    return upper_bound, lower_bound


# Normalization
def min_max_normalize(df):
    return (df - df.min()) / (df.max() - df.min())


def mean_normalize(df):
    return (df - df.mean()) / df.std()


def l2norm(x, pad=0.001):
    """
    Log2 normalization function

    :param float x: Input value
    :param int|float pad: Pad value (to handle zeros)
    :return: log2(x+1) normalized value
    :rtype: float
    """
    return np.log2(x + pad)
