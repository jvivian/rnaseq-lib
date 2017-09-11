import pandas as pd
import holoviews as hv

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
