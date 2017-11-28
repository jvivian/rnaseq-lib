import holoviews as hv

from rnaseq_lib.diff_exp import log2fc
from rnaseq_lib.tissues import subset_by_dataset
import numpy as np

df_cols = ['id', 'tissue', 'dataset', 'tumor', 'type']


def gene_curves(df, tissue, gene):
    # Subset dataframe for gene and tissue
    df = df[df.tissue == tissue]
    df = df[df_cols + [gene]].sort_values(gene, ascending=False)

    # Logscale gene for calculations
    df[gene] = df[gene].apply(lambda x: np.log2(x+1))

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
        tumor_mean = tumor.iloc[:int(len(tumor) * perc_tumor) - 1][gene].apply(lambda x: 2**x - 1).median()
        gtex_mean = gtex[gene].apply(lambda x: 2**x - 1).median()
        l2fc = log2fc(tumor_mean, gtex_mean)

        # Store
        records.append((tissue, exp, l2fc, perc_tumor, perc_normal, len(gtex), len(tumor), 'GTEx'))

    # Define dimensions
    tissue_dim = hv.Dimension('tissue', label='Tissue')
    ptumor_dim = hv.Dimension('percent_tumor', label='% Tumor')
    pnormal_dim = hv.Dimension('percent_normal', label='percent')
    l2fc_dim = hv.Dimension('l2fc', label='log2FC')
    exp_dim = hv.Dimension('expression', label='log2(x+1)')

    # First plot - Percentage of Normal Samples
    c1 = hv.Curve(data=df, kdims=[ptumor_dim],
                  vdims=[pnormal_dim, tissue_dim], group='Percentage of Normal Samples',
                  extents=(None, 0, None, 1))

    s1 = hv.Scatter(data=df, kdims=[ptumor_dim],
                    vdims=[pnormal_dim, tissue_dim], group='Percentage of Normal Samples')

    # Second Plot - Expression
    c2 = hv.Curve(data=df, kdims=[ptumor_dim],
                  vdims=[exp_dim, tissue_dim], group='Gene Expression',
                  extents=(None, 0, None, 16))

    s2 = hv.Scatter(data=df, kdims=[ptumor_dim],
                    vdims=[exp_dim, tissue_dim], group='Gene Expression')

    # Third Plot - Log2 Fold Change
    c3 = hv.Curve(data=df, kdims=[ptumor_dim],
                  vdims=[l2fc_dim, tissue_dim], group='Log2 Fold Change',
                  extents=(None, -0.5, None, 8))

    s3 = hv.Scatter(data=df, kdims=[ptumor_dim],
                    vdims=[l2fc_dim, tissue_dim], group='Log2 Fold Change')

    return (c1 * s1 + c2 * s2 + c3 * s3).cols(1)
