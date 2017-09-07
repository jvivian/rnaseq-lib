import numpy as np
import pandas as pd
from sklearn.manifold import TSNE

from tete import tete


def run_tete(df, num_dims=2, num_neighbors=20):
    """
    Runs t-ETE dimensionality reduction

    :param pd.DataFrame df: Dataframe or numpy array. Features need to be columns.
    :param int num_dims: Number of dimensions to reduce the array down to.
    :param int num_neighbors: Number of neighbors to use during iteration.
    :return: Reduced matrix
    :rtype: np.array
    """
    return tete(np.array(df), num_dims=num_dims, num_neighbs=num_neighbors)


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
