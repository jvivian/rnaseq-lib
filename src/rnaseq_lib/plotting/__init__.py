import pandas as pd
import holoviews as hv


def plot_boxplot(df,
                 samples_a, samples_b,
                 label_a, label_b,
                 feature,
                 norm_func=None,
                 value_label='counts', group_label='dataset'):
    """
    Return holoviews boxplot object for a samples by feature DataFrame

    :param pd.DataFrame df: Input DataFrame
    :param list[str] samples_a: list of samples for group A
    :param list[str] samples_b: list of samples for group B
    :param str label_a: Label for group A
    :param str label_b: Label for group B
    :param str feature: Feature (column label) to use
    :param func norm_func: Normalization function for dataframe
    :param str value_label: Label to use for values in boxplot
    :param str group_label: Label to use for groups in dataset
    :return: Holoviews boxplot object
    :rtype: hv.BoxWhisker
    """
    # Apply normalization function if provided
    if norm_func:
        df = df.apply(norm_func)

    plot = {value_label: df.loc[samples_a][feature].tolist() + df.loc[samples_b][feature].tolist(),
            group_label: [label_a for _ in samples_a] + [label_b for _ in samples_b]}

    return hv.BoxWhisker(pd.DataFrame.from_dict(plot), kdims=['dataset'], vdims=['counts'])
