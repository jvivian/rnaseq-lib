from __future__ import division

import holoviews as hv
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from rnaseq_lib.diff_exp import log2fc
from rnaseq_lib.dim_red import run_tsne, run_tete
from rnaseq_lib.plot.opts import *
from rnaseq_lib.tissues import subset_by_dataset


class Holoview:
    """
    Object for Holoviews plots of gene expression data. Created for use with Holomap and DynamicMap which cannot
    accept dataframes as arguments. This class circumvents that limitation by referencing
    the dataframe internally.
    """

    def __init__(self, df):
        """
        :param pd.DataFrame df: Dataframe containing metadata / expression values (Synapse.org: syn11515015)
        """
        self.df = df
        self.df_cols = ['tissue', 'tumor', 'type', 'labels']

        # Style attributes - used in conjunction with '.opts()'
        self._gene_curves_opts = gene_curves_opts
        self._gene_kde_opts = gene_kde_opts
        self._gene_distribution_opts = gene_distribution_opts
        self._gene_de_opts = gene_de_opts
        self._sample_count_opts = sample_count_opts
        self._l2fc_by_perc_samples_opts = l2fc_by_perc_samples_opts
        self._gene_de_heatmap_opts = gene_de_heatmap_opts
        self._de_concordance_opts = de_concordance_opts
        self._dist_with_iqr_bounds_opts = dist_with_iqr_bounds_opts
        self._dr_opts = dr_opts
        self._tissue_de_opts = tissue_de_opts
        # Hacky, but 5S_rRNA is the first gene in the annotation set
        self.samples = self.df.index.tolist()
        try:
            self._gene_start = self.df.columns.tolist().index('5S_rRNA')
            self.genes = self.df.columns[self._gene_start:].tolist()
            self.meta_cols = self.df.columns[:self._gene_start]
            self.samples = self.df.index.tolist()
            self.tissues = sorted(self.df.tissue.unique().tolist())
        except (ValueError, AttributeError) as e:
            print 'Missing attributes: \n{}'.format(e.message)
            self._gene_start = None
            self.genes = None
            self.meta_cols = None
            self.samples = None
            self.tissues = None

    # Internal methods
    def _subset(self, genes=None, tissue_subset=None):
        # Subset dataframe by gene
        df = self.df[self.df_cols + genes] if genes else self.df

        # Subset by tissues
        if tissue_subset:
            df = df[df.tissue.isin(tissue_subset)]
        return df

    def _gene_cutoff(self, gene, tissue, percent):
        # Subset dataframe by tissue and gene
        df = self._subset([gene], [tissue])

        # Subset by dataset
        tumor, normal, gtex = subset_by_dataset(df)

        # Calculate gene expression cutoffs for each dataset
        cutoffs = [x[gene].apply(self.l2norm).sort_values(ascending=False).iloc[int(len(x) * percent) - 1]
                   for x in [tumor, normal, gtex]]

        # Return mapping of dataset to cutoff
        return {x: y for x, y in zip(['tumor', 'normal', 'gtex'], cutoffs)}

    def _sample_counts_df(self, include_gtex=True):
        """
        Compute sample counts and returns dataframe

        :param bool include_gtex: If True, includes GTEx in returned sample counts
        :return: Sample counts for tissues and datasets
        :rtype: pd.DataFrame
        """
        # Subset by dataset
        tumor, normal, gtex = subset_by_dataset(self.df)

        dfs = [tumor, normal, gtex] if include_gtex else [tumor, normal]
        labels = ['tcga-tumor', 'tcga-normal', 'gtex'] if include_gtex else ['tcga-tumor', 'tcga-normal']

        # Count samples for each tissue
        records = []
        for tissue in sorted(self.df.tissue.unique()):
            for df, label in zip(dfs, labels):
                count = len(df[df.tissue == tissue])
                records.append([tissue, label, count])

        # Create dataframe
        return pd.DataFrame.from_records(records, columns=['Tissue', 'Label', 'Count']).sort_values(['Tissue', 'Label'])

    # Convenience methods
    @staticmethod
    def l2norm(x):
        """
        Log2 noramlization function for gene counts

        :param float x: Input value
        :return: log2(x+1) normalized value
        :rtype: float
        """
        return np.log2(x + 1)

    @staticmethod
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

    def perc_tumor_overexpressed(self, gene, tissue_subset=None):
        """
        Calculate the percent of tumor samples that overexpress a gene relative to the combined normal
        distribution of all samples (in tissue_subset) as well as compared to normals in the same tissue

        :param str gene: Gene (ex: ERBB2) to select
        :param list tissue_subset: List of tissues to subset by
        :return: Table of tissues and corresponding percentages pertaining to upper bound cutoff
        :rtype: pd.DataFrame
        """

        # Subset by gene and tissue
        df = self._subset(gene, tissue_subset)

        # Calculate upper and lower bounds across all tissues
        upper, lower = self.iqr_bounds(df[df.tumor == 'no'][gene].apply(self.l2norm))

        records = []
        for tissue in sorted(df.tissue.unique()):
            # Calclulate upper/lower bound for normal
            normals = df[(df.tumor == 'no') & (df.tissue == tissue)][gene].apply(self.l2norm)
            n_upper, n_lower = self.iqr_bounds(normals)

            # Calculate expression for tumor
            exp = df[(df.tissue == tissue) & (df.tumor == 'yes')][gene].apply(self.l2norm)

            # Calculate percentage cut offs
            perc = len([x for x in exp if x > upper]) / len(exp)
            n_perc = len([x for x in exp if x > n_upper]) / len(exp)
            records.append([tissue, perc, n_perc])

        # Save, sort, and return output dataframe
        df = pd.DataFrame.from_records(records, columns=['Tissue', 'Upper', 'T_Upper'])
        return df.sort_values('T_Upper', ascending=False)

    # Gene Plots
    def gene_kde(self, gene, tissue_subset=None, tumor=True, normal=True, gtex=True):
        """
        Returns KDE of gene expression (log2) for given tissue

        :param str gene: Gene (ex: ERBB2) to select
        :param list tissue_subset: List of tissues to subset by
        :param bool tumor: If True, include tumor samples
        :param bool normal: If True, include normal samples
        :param bool gtex: If True, include gtex samples
        :return: Returns holoviews Overlay object of gene KDE
        :rtype: hv.Overlay
        """
        # Subset dataframe by tissue and gene
        df = self._subset([gene], tissue_subset)

        # Subset by dataset
        t, n, g = subset_by_dataset(df)

        # Define x dimension for labeling
        x = hv.Dimension('Gene Expression', unit='log2(x+1)')

        # Create KDE objects for each tissue and dataset
        dists = []
        for tissue in df.tissue.unique():
            for label, dataset, flag in zip(['Tumor', 'GTEx', 'Normal'], [t, g, n], [tumor, gtex, normal]):
                if flag:
                    dists.append(hv.Distribution(dataset[dataset.tissue == tissue][gene].apply(self.l2norm),
                                                 kdims=[x], label='{}-{}'.format(label, tissue)))

        # Combine into Overlay object
        return hv.Overlay(dists, label='{} Expression'.format(gene)).opts(self._gene_kde_opts)

    def gene_distribution(self, gene, tissue_subset=None, types=None):
        """
        Box and Whisker expression distribution across tissues

        :param str gene: Gene (ex: ERBB2) to select
        :param list tissue_subset: List of tissues to subset by
        :param bool types: If True, uses tissue/cancer subtype instead of Tissue label
        :return: Returns holoviews BoxWhisker object
        :rtype: hv.BoxWhisker
        """
        # Subset dataframe by gene
        df = self._subset([gene], tissue_subset)

        # Normalize gene expression
        norm_exp = df[gene].apply(self.l2norm)

        # Subgroup labeling
        subgroup = df['type'] if types else df.tissue

        # Return grouped box and whiskers:
        return hv.BoxWhisker((subgroup, df['labels'], norm_exp), kdims=['Tissue', 'Dataset'],
                             vdims=[hv.Dimension('Gene Expression', unit='log2(x+1)')],
                             label='{} Expression'.format(gene)).opts(self._gene_distribution_opts)

    # Differential Expression
    def tissue_de(self, tissue, extents=None, tcga_normal=None, gene_labels=None):
        """
        Differential expression for a given tissue

        :param str tissue: Tissue to subset by
        :param tuple extents: xmin/ymin/xmax/ymax values
        :param bool tgca_normal: If True, uses TCGA normal for DE comparison
        :return: Scatterplot of DE
        :rtype: hv.Scatter
        """
        # Subset by tissue
        df = self._subset(genes=None, tissue_subset=[tissue])

        # Subset by dataset
        tumor, normal, gtex = subset_by_dataset(df)

        # Subset by genes
        t_genes = tumor[tumor.columns[self._gene_start:]]

        if tcga_normal:
            n_genes = normal[normal.columns[self._gene_start:]]
            label = 'Normal'
        else:
            n_genes = gtex[gtex.columns[self._gene_start:]]
            label = 'GTEx'

        # Calculate total expression
        exp = pd.concat([t_genes, n_genes]).apply(self.l2norm).median()

        # Calculate L2FC
        l2fc = log2fc(t_genes.median(), n_genes.median())

        # Define dimensions
        kdims = [hv.Dimension('exp', label='Gene Expression', unit='log2(x+1)')]
        vdim = hv.Dimension('l2fc', label='Log2 Fold Change', unit='log2(Tumor/{})'.format(label))

        plot = pd.DataFrame()
        plot['exp'] = exp
        plot['l2fc'] = l2fc
        plot['gene'] = exp.index
        plot.index = exp.index

        # Apply label column
        if gene_labels:
            label_vector = ['Other' for _ in plot.index]
            size_vector = [1 for _ in plot.index]
            for k, v in gene_labels.iteritems():
                for i in xrange(len(plot.index)):
                    if plot.index[i] in v:
                        label_vector[i] = k
                        size_vector[i] = 5
            plot['label'] = label_vector
            plot['size'] = size_vector
            vdims = [vdim] + ['gene', 'label', 'size']
        else:
            vdims = [vdim, 'gene']

        if extents:
            return hv.Scatter(plot, kdims=kdims, vdims=vdims, extents=extents).opts(self._tissue_de_opts)
        else:
            return hv.Scatter(plot, kdims=kdims, vdims=vdims).opts(self._tissue_de_opts)

    def gene_de(self, gene, tissue_subset=None, extents=None, tcga_normal=False):
        """
        Scatter plot of differential expression across all tissues

        :param str gene: Gene (ex: ERBB2) to select
        :param tuple extents: xmin/ymin/xmax/ymax values
        :param list tissue_subset: List of tissues to subset by
        :param bool tcga_normal: If True, will use TCGA normal for differential expression comparison
        :return: Scatterplot of values
        :rtype: hv.Scatter
        """
        # Subset dataframe by gene and tissue subset
        df = self._subset([gene], tissue_subset)

        # For each tissue, calculate L2FC and mean expression
        records = []
        for tissue in sorted(df.tissue.unique()):
            # Subset by dataset
            tumor, normal, gtex = subset_by_dataset(df[df.tissue == tissue])

            # Calculate tumor and normal expression and L2FC
            if tcga_normal:
                l2fc = log2fc(tumor[gene].median(), normal[gene].median())
                exp = pd.concat([tumor[gene], normal[gene]], axis=0).apply(self.l2norm).median()
                unit = 'log2(Tumor/Normal)'
            else:
                l2fc = log2fc(tumor[gene].median(), gtex[gene].median())
                exp = pd.concat([tumor[gene], gtex[gene]], axis=0).apply(self.l2norm).median()
                unit = 'log2(Tumor/GTEx)'

            # Store as record
            records.append((exp, l2fc, tissue))

        # Define dimensions of plot
        kdims = [hv.Dimension('Expression', label='Gene Expression', unit='log2(x+1)')]
        vdims = [hv.Dimension('L2FC', label='Fold Change', unit=unit), 'Tissue']

        # Create dataframe
        plot = pd.DataFrame.from_records(records, columns=['Expression', 'L2FC', 'Tissue'])

        if extents:
            return hv.Scatter(plot, kdims=kdims, vdims=vdims, extents=extents).opts(self._gene_de_opts)
        else:
            return hv.Scatter(plot, kdims=kdims, vdims=vdims).opts(self._gene_de_opts)

    def gene_de_kde(self, gene, tissue_subset=None, tcga_normal=False):
        """
        KDE of L2FC values for the tumor as compared to the normal

        :param str gene: Gene (ex: ERBB2) to select
        :param list tissue_subset: List of tissues to subset by
        :param bool tcga_normal: If True, use TCGA normal to for DE calc, otherwise use GTEx
        :return: Collection of Distribution objects
        :rtype: hv.Overlay
        """
        # Subset dataframe by gene and tissue subset
        df = self._subset([gene], tissue_subset)

        # Subset by dataset
        tumor, normal, gtex = subset_by_dataset(df)

        # Create X dimension
        xdim = hv.Dimension('Log2 Fold Change', unit='log2(a+1)/log2(b+1)')

        dists = []
        for tissue in df.tissue.unique():
            # Calculate mean expression for normal
            if tcga_normal:
                n = normal[normal.tissue == tissue][gene].median()
                label = 'Tumor-Normal-{}'.format(tissue)
            else:
                n = gtex[gtex.tissue == tissue][gene].median()
                label = 'Tumor-GTEx-{}'.format(tissue)

            # Calculate l2fc for each tumor sample and save
            l2fcs = []
            for i, row in tumor[tumor.tissue == tissue].iterrows():
                l2fcs.append(log2fc(row[gene], n))

            # Create distribution
            dists.append(hv.Distribution(l2fcs, kdims=[xdim], label=label))

        return hv.Overlay(dists, label='{} Expression'.format(gene)).opts(self._gene_kde_opts)

    def l2fc_by_perc_samples(self, gene, tissue_subset=None, tcga_normal=False, l2fc_cutoff=2):
        """
        Calculate the percentage of samples greater than a range of log2 fold change values

        :param str gene: Gene (ex: ERBB2) to select
        :param list tissue_subset: List of tissues to subset by
        :param bool tcga_normal: If True, use TCGA normal to for DE calc, otherwise use GTEx
        :param float l2fc_cutoff: Specifies the L2FC cutoff to draw a Spike object
        :return: Collection of Curve objects
        :rtype: hv.Overlay
        """
        # Subset dataframe by gene and tissue subset
        df = self._subset([gene], tissue_subset)

        # Subset by dataset
        tumor, normal, gtex = subset_by_dataset(df)

        # Create X dimension
        xdim = hv.Dimension('Log2 Fold Change', unit='log2(a+1)/log2(b+1)')
        ydim = hv.Dimension('Tumor Samples With Greater L2FC', unit='%')

        # Calculate % samples over a given l2fc
        curves = []
        label = ''
        for tissue in sorted(df.tissue.unique()):
            # Calculate mean expression for normal
            if tcga_normal:
                n = normal[normal.tissue == tissue][gene].median()
                label = 'Tumor-Normal'
            else:
                n = gtex[gtex.tissue == tissue][gene].median()
                label = 'Tumor-GTEx'

            # Calculate l2fc for each tumor sample and save
            l2fcs = []
            for i, row in tumor[tumor.tissue == tissue].iterrows():
                l2fcs.append(log2fc(row[gene], n))

            # Calculate percentage samples over l2fc
            percentages = {}
            l2fc_range = [x * 0.1 for x in xrange(0, int(np.ceil(max(l2fcs) * 10)))]
            for l2fc in l2fc_range:
                percentages[l2fc] = len([x for x in l2fcs if x >= l2fc]) / len(l2fcs) * 100

            # Create line object
            curves.append(hv.Area(percentages, kdims=[xdim], vdims=[ydim], label=tissue))

        # Return curves along with a Spikes object at the l2fc cutoff
        overlay = hv.Overlay(curves + [hv.Spikes([l2fc_cutoff])], label='{} {} Expression'.format(label, gene))
        return overlay.opts(self._l2fc_by_perc_samples_opts)

    def gene_de_heatmap(self, genes, tissue_subset=None, tcga_normal=False):
        """
        Heatmap of gene log2 fold change

        :param str gene: Gene (ex: ERBB2) to select
        :param list tissue_subset: List of tissues to subset by
        :param bool tcga_normal: If True, use TCGA normal to for DE calc, otherwise use GTEx
        :return: DE Heatmap of genes for tissue subset
        :rtype: hv.HeatMap
        """
        # Subset dataframe by genes
        df = self.df[self.df_cols + genes]

        # Subset by tissues
        if tissue_subset:
            df = df[df.tissue.isin(tissue_subset)]

        # Subset by dataset
        tumor, normal, gtex = subset_by_dataset(df)

        # For each tissue/gene, calculate L2FC
        records = []
        for tissue in sorted(df.tissue.unique()):
            for gene in genes:
                # Calculate mean expression for normal
                if tcga_normal:
                    n = normal[normal.tissue == tissue][gene].median()
                else:
                    n = gtex[gtex.tissue == tissue][gene].median()

                # Calculate expression for tumor and compute l2fc
                t = tumor[tumor.tissue == tissue][gene].median()
                l2fc = log2fc(t, n)

                records.append([tissue, gene, l2fc])

        # Create dataframe and define dimensions
        df = pd.DataFrame.from_records(records, columns=['Tissue', 'Gene', 'L2FC']).sort_values('Tissue')

        return hv.HeatMap(df, kdims=['Gene', 'Tissue'], vdims=['L2FC']).opts(self._gene_de_heatmap_opts)

    def tissue_top_de_genes(self, tissue):
        # Create DE objects to get data
        gtex = self.tissue_de(tissue).data
        tcga = self.tissue_de(tissue, tcga_normal=True).data

        intervals = [10, 100, 500, 1000, 5000, 10000, len(self.genes)]

        # Calculate maximum arange for plot
        reg_line_arange = gtex[gtex.exp > gtex.exp.median()].sort_values('l2fc', ascending=False).l2fc.tolist()

        # Top DE genes with high expression
        hmaps = {}
        for i in intervals:
            x = gtex[gtex.exp > gtex.exp.median()].sort_values('l2fc', ascending=False).l2fc.tolist()[:i]
            y = tcga[tcga.exp > tcga.exp.median()].sort_values('l2fc', ascending=False).l2fc.tolist()[:i]

            scatter = hv.Scatter((x, y), kdims=['GTEx L2FC'], vdims=['TCGA L2FC'])
            reg_line = hv.Curve(self.regression_line(x, y, arange=reg_line_arange))
            pearson_r = round(pearsonr(x, y)[0], 2)

            title = 'R: {}'.format(pearson_r)
            hmaps[i] = hv.Overlay([scatter, reg_line]).relabel(title)
        return hv.HoloMap(hmaps, kdims='Num_Genes').relabel('Top DE Gene L2FC in {}'.format(tissue))

    # Misc plots
    def dist_with_iqr_bounds(self, ys, kdim):
        """
        Creates distribution object with IQR bounds

        :param list ys: List of values to calculate IQR and bounds
        :param str kdims: K-dimension label for distribution
        :return: Distribution with IQR bounds
        :rtype: hv.Overlay
        """
        # Calculate IQR and outlier bounds
        q25, q75 = np.percentile(ys, [25, 75])
        upper, lower = self.iqr_bounds(ys)

        # Return dist with spikes
        return hv.Overlay([hv.Distribution(ys, kdims=[kdim]),
                           hv.Spikes([q25, q75]),
                           hv.Spikes([lower, upper])]).opts(self._dist_with_iqr_bounds_opts)

    @staticmethod
    def regression_line(x, y, arange=None):
        """
        Returns x/y vectors of a regression line for 2D input

        :param np.array x: Vector of x values
        :param np.array y: Vector of y values
        :param np.array arange: Provide a custom arange to generate regression line
        :return: Regression line vectors
        :rtype: tuple(np.array, np.array)
        """
        m, b = np.polyfit(x, y, 1)
        reg_x = np.arange(min(arange), max(arange)) if arange else np.arange(min(x), max(x))
        return reg_x, m * reg_x + b

    @staticmethod
    def path_box(xmin, xmax, ymin, ymax, color=None):
        """
        Returns rectangular Path object for a given set of x/y coordinates

        :param float xmin: xmin of box
        :param float xmax: xmax of box
        :param float ymin: ymin of box
        :param float ymax: ymax of box
        :param str color: Set the color of the Path object
        :return: Rectangular path object
        :rtype: hv.Path
        """
        path = [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax), (xmin, ymin)]
        if color:
            return hv.Path([path]).opts(dict(Path=dict(style=dict(color=color))))
        else:
            return hv.Path([path])

    def highlight_points(self, xs, ys, size=0.1, color=None, hidden_buffer_box=True):
        """
        Returns a rectangular Path object for a set of points

        :param list|float xs: List of x coordinates or a single x coord
        :param list|float ys: List of y coordinates or a single y coord
        :param float size: Margin around xmin,xmax,ymin,ymax of points
        :param str color: Set the color of the Path object
        :param bool hidden_buffer_box: Adds a transparent larger frame around the Path object to improve plot margins
        :return: Rectangular Path object
        :rtype: hv.Path
        """
        # If a set of single points
        if isinstance(xs, (int, float)) and isinstance(ys, (int, float)):
            xs, ys = [xs], [ys]

        # Collect mins nad maxes from all points
        xmin, xmax, ymin, ymax = min(xs), max(xs), min(ys), max(ys)

        # Add margins
        xmin, xmax, ymin, ymax = xmin - size, xmax + size, ymin - size, ymax + size

        # Create Path object
        plot = self.path_box(xmin, xmax, ymin, ymax, color=color)

        # If hidden_buffer_box is enabled
        if hidden_buffer_box:
            xmin, xmax, ymin, ymax = xmin - size, xmax + size, ymin - size, ymax + size
            hbb = self.path_box(xmin, xmax, ymin, ymax).opts(dict(Path=dict(style=dict(alpha=0))))
            return plot * hbb
        else:
            return plot

    def gene_curves(self, gene, tissue):
        """
        Returns set of 3 plots for tissue / gene given a dataframe of metadata and expression values

        :param str gene: Gene (ex: ERBB2) to select
        :param str tissue: Tissue (ex: Breast) to select
        :return: Returns holoviews Layout object containing 3 plots for selected Tisssue / Gene
        :rtype: hv.Layout
        """
        # Subset dataframe for gene and tissue
        df = self._subset([gene], [tissue])

        # Logscale gene for calculations
        df[gene] = df[gene].apply(self.l2norm)

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

    def sample_counts(self, include_gtex=True):
        """
        Bargraph of tissues grouped by dataset

        :param bool include_gtex: If True, GTEx is included in sample counts
        :return: Bargraph of sample counts
        :rtype: hv.Bars
        """
        df = self._sample_counts_df(include_gtex=include_gtex).sort_values(['Tissue', 'Label'], ascending=[True, False])

        # Return Bars object of sample counts
        return hv.Bars(df, kdims=['Tissue', 'Label'], vdims=['Count'],
                       label='Sample Counts for TCGA and GTEx').opts(self._sample_count_opts)

    def differential_expression_tissue_concordance(self, tissue_subset=None, tcga=True, gtex=True):
        """
        Categorical scatterplot of concordance between tissues for gene differential expression

        :param list tissue_subset: List of tissues to subset by
        :return: Heatmap of differential expression comparison across tissue
        :rtype: hv.HeatMap
        """

        df = self._subset(genes=None, tissue_subset=tissue_subset)
        records = []
        for tissue1 in sorted(df.tissue.unique()):

            # Subset by tissue then break apart by dataset
            t, g, n = subset_by_dataset(df[df.tissue == tissue1])

            # If there are both normal and gtex samples
            if len(g) > 0 and len(n) > 0:

                # Calculate gene expression average for tumor samples
                i = df.columns.tolist().index('OR4F5')  # Hacky, but OR4F5 is the first gene in the dataframe
                tmed = t[t.columns[i:]].median()
                nmed = n[n.columns[i:]].median()
                master_tn = log2fc(tmed, nmed)

                # Iterate over all other tissues to get PearsonR of L2FC
                for tissue2 in sorted(df.tissue.unique()):

                    # Subset by second tissue
                    _, g, n = subset_by_dataset(df[df.tissue == tissue2])

                    # If there are GTEx samples, calculate PearsonR to master_tn
                    if gtex and len(g) > 0:
                        gmed = g[g.columns[i:]].median()
                        tg = log2fc(tmed, gmed)
                        records.append((tissue1, '{}-GTEx'.format(tissue2), round(pearsonr(master_tn, tg)[0], 2)))

                    # If there are normal samples, calculate PearsonR to master_tn
                    if tcga and len(n) > 0:
                        nmed = n[n.columns[i:]].median()
                        tn = log2fc(tmed, nmed)
                        records.append((tissue1, '{}-TCGA'.format(tissue2), round(pearsonr(master_tn, tn)[0], 2)))


        # Construct dataframe
        df = pd.DataFrame.from_records(records, columns=['Tissue-Tumor/Normal', 'Tissue-Normal', 'PearsonR'])

        # Return HeatMap object
        return hv.HeatMap(df, kdims=['Tissue-Tumor/Normal', 'Tissue-Normal'], vdims=['PearsonR'],
                          label='Differential Expression Gene Concordance (PearsonR)').opts(self._de_concordance_opts)

    # Dimensionality Reduction
    def trimap(self, genes, tissue_subset=None, num_neighbors=50):
        """
        Dimensionality reduction via Trimap

        :param list(str) genes: List of genes to subset by
        :param list(str) tissue_subset: List of tissues to subset by
        :param int num_neighbors: Hyperparameter for trimap
        :return: Scatterplot of dimensionality reduction
        :rtype: hv.Scatter
        """
        # Subset dataframe by genes (keeping some metadata)
        df = self.df[self.df_cols + genes].sort_values('tissue')

        # Subset by tissues
        if tissue_subset:
            df = df[df.tissue.isin(tissue_subset)]

        # Run Trimap (used to be called t-ETE)
        z = run_tete(df[genes], num_dims=2, num_neighbors=num_neighbors)

        # Add results to dataframe
        df['x'] = z[:, 0]
        df['y'] = z[:, 1]

        return hv.Scatter(df, kdims=['x'], vdims=['y'] + self.df_cols).opts(self._dr_opts)

    def tsne(self, genes, tissue_subset=None, perplexity=50, learning_rate=1000):
        """
        Dimensionality reduction via t-SNE

        :param list(str) genes: List of genes to subset by
        :param list(str) tissue_subset: List of tissues to subset by
        :param int perplexity: Hyperparamter for t-SNE
        :param int learning_rate: Hyperparamter for t-SNE
        :return: Scatterplot of dimensionality reduction
        :rtype: hv.Scatter
        """
        # Subset dataframe by genes (keeping some metadata)
        df = self.df[self.df_cols + genes].sort_values('tissue')

        # Subset by tissues
        if tissue_subset:
            df = df[df.tissue.isin(tissue_subset)]

        # Run t-SNE
        z = run_tsne(df[genes], num_dims=2, perplexity=perplexity, learning_rate=learning_rate)

        # Add results to dataframe
        df['x'] = z[:, 0]
        df['y'] = z[:, 1]

        return hv.Scatter(df, kdims=['x'], vdims=['y'] + self.df_cols).opts(self._dr_opts)


def disable_logo(plot, element):
    plot.state.toolbar.logo = None

#hv.plotting.bokeh.ElementPlot.finalize_hooks.append(disable_logo)