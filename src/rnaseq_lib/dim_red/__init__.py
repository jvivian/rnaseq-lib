import numpy as np
import pandas as pd
from sklearn.manifold import TSNE

from trimap import trimap


def run_trimap(df, num_dims=2, kin=50, kout=5, krand=5, eta=10000.0):
    """
    Runs t-ETE dimensionality reduction

    :param pd.DataFrame df: Dataframe or numpy array. Features need to be columns.
    :param int num_dims: Number of dimensions to reduce the array down to.
    :param int kin: Number of k-Nearest Neighbor points
    :param int kout: Number of outliers (num triplets per point = kin * kout)
    :param int krand: Number of random triplets per point
    :param float eta: Initial learning rate
    :return: Reduced matrix
    :rtype: np.array
    """
    return trimap(np.array(df), num_dims=num_dims, kin=kin, kout=kout, krand=krand, eta=eta)


def run_tsne(df, num_dims, perplexity=30, learning_rate=200):
    """
    Runs t-SNE dimensionality reduction

    :param pd.DataFrame df: Dataframe or numpy array. Features need to be columns.
    :param int num_dims: Number of dimensions to reduce the array down to.
    :param int perplexity: Perplexity hyperparameter for t-SNE
    :param int learning_rate: Learning Rate hyperparameter for t-SNE
    :return: Reduced matrix
    :rtype: np.array
    """
    return TSNE(n_components=num_dims,
                perplexity=perplexity,
                learning_rate=learning_rate).fit_transform(np.array(df))
