import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import norm
from sklearn import mixture

from rnaseq_lib.math.dists import name_from_dist, DISTRIBUTIONS


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


def softmax(df):
    """
    Normalizes columns to sum to 1

    :param pd.DataFrame df: Dataframe to normalize
    :return: Normalized DataFrame
    :rtype: pd.DataFrame
    """
    return df.divide(df.sum())


def l2norm(x, pad=0.001):
    """
    Log2 normalization function

    :param float x: Input value
    :param int|float pad: Pad value (to handle zeros)
    :return: log2(x+1) normalized value
    :rtype: float
    """
    return np.log2(x + pad)


# Distributions
def run_ks(source_dist, dists=DISTRIBUTIONS):
    """
    Runs Kolmogorov-Smirnov test for the provided source distribution against provided scipy distribution funcs

    :param np.array source_dist: Distribution to test
    :param list(func) dists: List of scipy.stats distributions to test. Defaults to list containing most.
    :return: Dataframe of KS-test results
    :rtype: pd.DataFrame
    """
    rows = []
    for dist in dists:
        kstat, pval = stats.kstest(source_dist, name_from_dist(dist), args=dist.fit(source_dist))
        rows.append((name_from_dist(dist), kstat, pval))
    return pd.DataFrame(rows, columns=['Name', 'KS-stat', 'Pvalue']).sort_values('KS-stat')


def find_gaussian_intersection(m1, m2, std1, std2):
    """
    Given parameters for two gaussian distributions, identify the intersection(s)

    :param float m1: Mean for first Gaussian
    :param float m2: Mean for second Gaussian
    :param float std1: Standard deviation for first Gaussian
    :param float std2: Standard deviation for second Gaussian
    :return: Intersection(s) between Gaussian distributions
    :rtype: list(float,)
    """
    # Define systems of equations
    m1, m2, std1, std2 = float(m1), float(m2), float(std1), float(std2)
    a = 1.0 / (2 * std1 ** 2) - 1.0 / (2 * std2 ** 2)
    b = m2 / (std2 ** 2) - m1 / (std1 ** 2)
    c = m1 ** 2 / (2 * std1 ** 2) - m2 ** 2 / (2 * std2 ** 2) - np.log(std2 / std1)

    # Return intersection between means
    mean_min, mean_max = sorted([m1, m2])

    # Only return the intersection if one exists between the means
    roots = [round(x, 2) for x in np.roots([a, b, c])]
    inter = [x for x in np.roots([a, b, c]) if mean_min < x < mean_max]
    if len(inter) == 0:
        return roots
    else:
        return inter


def overlay_gmm_to_hist(source_dist, figsize=(12, 4), color='red'):
    """
    Given a source distribution, fit a 2-component Gaussian mixture model and return plot

    :param np.array source_dist: Source distribution
    :param tuple(int, int) figsize: Figure size
    :return: Fig object and cutoff value
    :rtype: float
    """
    # Fit GMM
    gmm = mixture.GaussianMixture(n_components=2).fit(pd.DataFrame(source_dist))
    m1, m2 = gmm.means_
    std1, std2 = gmm.covariances_

    # Identify intersection between the two Gaussians
    cutoffs = find_gaussian_intersection(m1, m2, std1, std2)

    # Plot source data
    plt.subplots(figsize=figsize)
    plt.hist(source_dist, density=True, alpha=0.25, bins=50, label='Tumor', color=color)

    # Plot Gaussian fits and intersection
    x = np.linspace(min(source_dist), max(source_dist), len(source_dist))
    plt.plot(x, *norm.pdf(x, m1, std1), label='u={}, o={}'.format(round(m1, 1), round(std1, 1)))
    plt.plot(x, *norm.pdf(x, m2, std2), label='u={}, o={}'.format(round(m2, 1), round(std2, 1)))
    # Add intersectional lines
    for cutoff in cutoffs:
        plt.vlines(cutoff, *plt.ylim(), label='Cutoff: {}'.format(cutoff), color='red', linestyles='--')
    plt.legend()
    if len(cutoffs) == 1:
        return cutoffs[0]
    else:
        return max(cutoffs)
