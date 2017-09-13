from collections import defaultdict

import holoviews as hv
import numpy as np
import pandas as pd

from rnaseq_lib.plotting import plot_boxplot
from rnaseq_lib.tissues import get_gtex_samples, return_samples
from rnaseq_lib.tissues import get_normal_samples
from rnaseq_lib.tissues import get_tumor_samples


def gene_expression_boxplot(df, tissue, gene, normal=False):
    """
    Returns a holoviews box and whisker object of a gene for a tissue

    :param pd.DataFrame df: Gene by Samples DataFrame
    :param str tissue: Tissue to use
    :param str gene: Gene to use
    :param bool normal: If True, will include TCGA normal samples
    :return: Boxplot of gene for tissue
    :rtype: hv.BoxWhisker
    """
    plot_info = {'GTEx': get_gtex_samples(tissue),
                 'Tumor': get_tumor_samples(tissue)}

    if normal:
        plot_info['Normal'] = get_normal_samples(tissue)

    return plot_boxplot(df, plot_info=plot_info, feature=gene, title=tissue, norm_func=lambda x: np.log2(x + 1))


def get_tcga_gtex_label_info():
    """
    :return: Returns TCGA and GTEx label, tissue, and samples information
    :rtype: defaultdict
    """
    samples = return_samples()
    info = defaultdict(list)
    for tissue in sorted(samples):
        info['label'].extend(['Tumor' if x.endswith('-01')
                              else 'Normal' if x.endswith('-11') else 'GTEx' for x in samples[tissue]])
        info['sample'].extend([x for x in samples[tissue]])
        info['tissue'].extend([tissue for _ in samples[tissue]])
    return info
