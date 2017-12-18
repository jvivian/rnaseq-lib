import holoviews as hv
import numpy as np
import pandas as pd

from rnaseq_lib.diff_exp import log2fc
from rnaseq_lib.dim_red import run_tsne, run_tete
from rnaseq_lib.plot.opts import gene_curves_opts, gene_kde_opts, gene_distribution_opts, gene_de_opts
from rnaseq_lib.tissues import subset_by_dataset


class Holoview:
    """
    Object for Holoviews plots. Created for use with Holomap and DynamicMap which cannot
    accept dataframes as arguments. This class circumvents that limitation by referencing
    the dataframe internally.
    """

    def __init__(self, df):
        """
        :param pd.DataFrame df: Dataframe containing metadata / expression values (Synapse.org: syn11515015)
        """
        self.df = df
        self.df_cols = ['id', 'tissue', 'dataset', 'tumor', 'type']

        # Style attributes - used in conjunction with '.opts()'
        self.gene_curves_opts = gene_curves_opts
        self.gene_kde_opts = gene_kde_opts
        self.gene_distribution_opts = gene_distribution_opts
        self.gene_de_opts = gene_de_opts

    @staticmethod
    def l2norm(x):
        """
        Log2 noramlization function for gene counts

        :param float x: Input value
        :return: log2(x+1) normalized value
        :rtype: float
        """
        return np.log2(x + 1)
        df = self.df[self.df_cols + [gene]].sort_values(gene, ascending=False)
        return df[df.tissue == tissue]

    def gene_kde(self, gene, tissue):
        """
        Returns KDE of gene expression (log2) for given tissue

        :param str gene: Gene (ex: ERBB2) to select
        :param str tissue: Tissue (ex: Breast) to select
        :return: Returns holoviews Overlay object of gene KDE
        :rtype: hv.Overlay
        """
        # Subset dataframe by tissue and gene
        df = self.subset(tissue, gene)

        # Subset by dataset
        tumor, normal, gtex = subset_by_dataset(df)

        # Define x dimension for labeling
        x = hv.Dimension('Gene Expression', unit='log2(x+1)')

        # Create KDE objects
        t = hv.Distribution(tumor[gene].apply(lambda x: np.log2(x + 1)), kdims=[x], label='Tumor-{}'.format(tissue))
        g = hv.Distribution(gtex[gene].apply(lambda x: np.log2(x + 1)), kdims=[x], label='GTEx-{}'.format(tissue))

        return hv.Overlay([t, g], label='{} Expression'.format(gene))

    def multiple_tissue_gene_kde(self, gene, *tissues):
        return hv.Overlay([self.gene_kde(gene, t) for t in tissues],
                          label='{} Expression'.format(gene))

    def gene_distribution(self, gene, extents=None):
        # Subset dataframe by gene
        df = self.df[self.df_cols + [gene]].sort_values(gene, ascending=False)

        # Normalize gene expression
        df[gene] = df[gene].apply(lambda x: np.log2(x + 1))

        # Subset for Tumor and GTEx
        df = df[((df.tumor == 'yes') | (df.dataset == 'gtex'))]

        # return grouped box and whiskers:
        return hv.BoxWhisker((df.tissue, df.dataset, df[gene]), kdims=['tissue', 'dataset'], vdims='gene')

    def gene_DE(self, gene, extents=None):
        # Subset dataframe by gene
        df = self.df[self.df_cols + [gene]].sort_values(gene, ascending=False)

        # Subset by dataset
        tumor, normal, gtex = subset_by_dataset(df)

        # For each tissue, calculate L2FC and mean expression
        records = []
        for tissue in sorted(df.tissue.unique()):
            exp = df[(df.tissue==tissue) & ((df.tumor == 'yes') | (df.dataset == 'gtex'))][gene].apply(lambda x: np.log2(x + 1)).median()

            t = tumor[tumor.tissue == tissue][gene].median()
            g = gtex[gtex.tissue == tissue][gene].median()
            l2fc = log2fc(t, g)

            records.append((exp, l2fc, tissue))

        kdims = ['Expression']
        vdims = ['L2FC', 'Tissue']

        plot = pd.DataFrame.from_records(records, columns=kdims + vdims)

        if extents:
            return hv.Scatter(plot, kdims=kdims, vdims=vdims, extents=extents)
        else:
            return hv.Scatter(plot, kdims=kdims, vdims=vdims)

    def gene_curves(self, gene, tissue):
        """
        Returns set of 3 plots for tissue / gene given a dataframe of metadata and expression values

        :param str gene: Gene (ex: ERBB2) to select
        :param str tissue: Tissue (ex: Breast) to select
        :return: Returns holoviews Layout object containing 3 plots for selected Tisssue / Gene
        :rtype: hv.Layout
        """
        # Subset dataframe for gene and tissue
        df = self.subset(tissue, gene)

        # Logscale gene for calculations
        df[gene] = df[gene].apply(lambda x: np.log2(x + 1))

        # Subset by dataset
        tumor, normal, gtex = subset_by_dataset(df)

        # Get values for plot
        records = []
        for perc_tumor in [x * 0.1 for x in xrange(1, 11)]:
            # Get log2 expression value for top x% tumor samples
            exp = float(tumor.iloc[int(len(tumor) * perc_tumor) - 1][gene])

            # Get percentage of samples in GTEx
            perc_normal = (len(gtex[gtex[gene] > exp]) * 1.0) / len(gtex)

            # Compute L2FC for tumor sample subset vs GTEx
            tumor_mean = tumor.iloc[:int(len(tumor) * perc_tumor) - 1][gene].apply(lambda x: 2 ** x - 1).median()
            gtex_mean = gtex[gene].apply(lambda x: 2 ** x - 1).median()
            l2fc = log2fc(tumor_mean, gtex_mean)

            # Store
            records.append((tissue, exp, l2fc, perc_tumor, perc_normal, len(gtex), len(tumor), 'GTEx'))

        # Create dataframe from records
        info = pd.DataFrame.from_records(records, columns=['tissue', 'expression',
                                                           'l2fc',
                                                           'percent_tumor',
                                                           'percent_normal',
                                                           'num_normals', 'num_tumors',
                                                           'normal_dataset'])

        # Define dimensions
        tissue_dim = hv.Dimension('tissue', label='Tissue')
        ptumor_dim = hv.Dimension('percent_tumor', label='% Tumor')
        pnormal_dim = hv.Dimension('percent_normal', label='percent')
        l2fc_dim = hv.Dimension('l2fc', label='log2FC')
        exp_dim = hv.Dimension('expression', label='log2(x+1)')

        # First plot - Percentage of Normal Samples
        c1 = hv.Curve(data=info, kdims=[ptumor_dim],
                      vdims=[pnormal_dim, tissue_dim], group='Percentage of Normal Samples',
                      extents=(None, 0, None, 1))

        s1 = hv.Scatter(data=info, kdims=[ptumor_dim],
                        vdims=[pnormal_dim, tissue_dim], group='Percentage of Normal Samples')

        # Second Plot - Expression
        c2 = hv.Curve(data=info, kdims=[ptumor_dim],
                      vdims=[exp_dim, tissue_dim], group='Gene Expression',
                      extents=(None, 0, None, 16))

        s2 = hv.Scatter(data=info, kdims=[ptumor_dim],
                        vdims=[exp_dim, tissue_dim], group='Gene Expression')

        # Third Plot - Log2 Fold Change
        c3 = hv.Curve(data=info, kdims=[ptumor_dim],
                      vdims=[l2fc_dim, tissue_dim], group='Log2 Fold Change',
                      extents=(None, -0.5, None, 8))

        s3 = hv.Scatter(data=info, kdims=[ptumor_dim],
                        vdims=[l2fc_dim, tissue_dim], group='Log2 Fold Change')

        return (c1 * s1 + c2 * s2 + c3 * s3).cols(1)


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
