import holoviews as hv
import numpy as np
import pandas as pd
from rnaseq_lib.dim_red import run_tsne, run_tete
from rnaseq_lib.utils import flatten


def plot_boxplot(df,
                 plot_info,
                 feature,
                 norm_func=None,
                 title=None,
                 value_label='counts', group_label='dataset'):
    """
    Return holoviews boxplot object for a "samples by feature" DataFrame

    :param pd.DataFrame df: Input DataFrame
    :param dict(str, list(str)) plot_info: Dict in the form "Label: [Samples]"
    :param str feature: Feature (column label) to use
    :param func norm_func: Normalization function for dataframe
    :param str title: Title of plot
    :param str value_label: Label to use for values in boxplot
    :param str group_label: Label to use for groups in dataset
    :return: Holoviews boxplot object
    :rtype: hv.BoxWhisker
    """
    # Apply normalization function if provided
    if norm_func:
        df = df.apply(norm_func)

    # Define group label
    group = []
    for label in sorted(plot_info):
        group.extend([label for _ in plot_info[label]])

    # Create dictionary with plot info
    plot = {value_label: flatten([df.loc[plot_info[x]][feature].tolist() for x in sorted(plot_info)]),
            group_label: group}

    # Return Holoviews BoxWhisker object
    return hv.BoxWhisker(pd.DataFrame.from_dict(plot), kdims=['dataset'], vdims=['counts'], group=title)


def tsne_of_dataset(df, title, perplexity=50, learning_rate=1000, plot_info=None):
    """
    t-SNE plot of a dataset

    :param pd.DataFrame df: Samples by features DataFrame
    :param str title: Title of plot
    :param int perplexity: Perplexity hyperparamter for t-SNE
    :param int learning_rate: Learning rate hyperparameter for t-SNE
    :param dict plot_info: Additional information to include in plot
    :return: Holoviews scatter object
    :rtype: hv.Scatter
    """
    z = run_tsne(df, num_dims=2, perplexity=perplexity, learning_rate=learning_rate)
    return _scatter_dataset(z=z, title=title, info=plot_info)


def tete_of_dataset(df, title, num_neighbors=30, plot_info=None):
    """
    t-ETE plot of a dataset

    :param pd.DataFrame df: Samples by features DataFrame
    :param str title: Title of plot
    :param int num_neighbors: Number of neighbors in t-ETE algorithm
    :param dict plot_info: Additional information to include in plot
    :return: Holoviews scatter object
    :rtype: hv.Scatter
    """
    z = run_tete(df, num_dims=2, num_neighbors=num_neighbors)
    return _scatter_dataset(z, title=title, info=plot_info)


def _scatter_dataset(z, title, info=None):
    """
    Internal function for scattering dataset

    :param np.array z: An [n x 2] matrix of values to plot
    :param dict info: Additional info for plotting. Lengths of values must match x and y vectors derived from z
    """
    # Collect information for plotting
    if info is None:
        info = dict()

    info['x'] = z[:, 0]
    info['y'] = z[:, 1]

    # Return Holoviews Scatter object
    return hv.Scatter(pd.DataFrame.from_dict(info),
                      kdims=['x'],
                      vdims=['y'] + [x for x in info.keys() if not x == 'x' and not x == 'y'],
                      group=title)


def plot_deseq2(df, info):
    z = np.array(df[['baseMean', 'log2FoldChange']])
    z[:, 0] = map(lambda x: np.log2(x + 1), z[:, 0])
    info['x'] = z[:, 0]
    info['y'] = z[:, 1]
    info['Name'] = df.index

    return hv.Scatter(pd.DataFrame.from_dict(info),
                      kdims=[('x', 'Log2 Mean')],
                      vdims=[('y', 'Log2 Fold Change')] + [k for k in info.keys() if k not in ['x', 'y']],
                      group='Lung')
