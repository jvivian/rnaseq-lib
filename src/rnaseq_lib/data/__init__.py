import os
import pickle

import pandas as pd

_cwd = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


def add_metadata_to_exp(exp, met):
    """
    Adds metadata to the expression dataframe and returns the combined object

    :param pd.DataFrame exp: Expression dataframe
    :param pd.DataFrame met: Metadata dataframe
    :return: Expression dataframe with added metadata columns
    :rtype: pd.DataFrame
    """
    # Copy genes from expression dataframe
    genes = exp.columns.tolist()

    # Remove duplicates from metadata
    samples = [x for x in exp.index if x in met.id]
    met = met[met.id.isin(samples)].drop_duplicates('id')

    # Ensure index dims are the same length
    assert len(exp) == len(met), 'Expression dataframe and metadata do not match index lengths'

    # Add metadata and return resorted dataframe
    df = exp.copy()
    df['id'] = met.id
    df['tissue'] = met.tissue
    df['type'] = met.type
    df['tumor'] = met.tumor
    df['label'] = label_vector_from_samples(df.index)
    return df[['id', 'tissue', 'type', 'label', 'tumor'] + genes]


def label_vector_from_samples(samples):
    """
    Produce a vector of labels for the sample vector provided

    :param list(str) samples: List of samples to derive labels for
    :return: Label vector
    :rtype: list(str)
    """
    vector = []
    for x in samples:
        if x.startswith('TCGA'):
            if x.endswith('11'):
                vector.append('tcga-normal')
            elif x.endswith('01'):
                vector.append('tcga-tumor')
            else:
                vector.append('tcga-other')
        else:
            vector.append('gtex')
    return vector


def map_genes(genes, strict=True):
    """
    Maps gene IDs to gene names

    :param list genes: ENSEMBL gene IDs to be mapped to gene names
    :param bool strict: If true, raies a KeyError if gene is not found in the gene_map
    :return: Mapped genes
    :rtype: list
    """
    gene_map = load_gene_map()
    if strict:
        return [gene_map[x.split('.')[0]] for x in genes]
    else:
        mapped = []
        for g in genes:
            try:
                mapped.append(gene_map[g.split('.')[0]])
            except KeyError:
                mapped.append(g)
        return mapped


def get_ucsf_subset(df, strict=False):
    """
    Subset UCSF dataframe and return.

    :param pd.DataFrame df: Input Dataframe in the format of "Genes by Samples"
    :param bool strict: If True, raises an error if gene is unmapped
    :return: Subset of Dataframe that only includes UCSF genes
    :rtype: pd.DataFrame
    """
    df.index = map_genes(df.index, strict=strict)

    ucsf_genes = load_ucsf_genes()
    ucsf_genes = [x for x in ucsf_genes if x in df.index]

    return df.loc[ucsf_genes]


def yield_sample_type(samples, kind='tumor'):
    """
    Yield samples based on what kind is being requested

    :param list(str) samples: List of samples
    :param str kind: What to return. 'tumor' and 'normal' for TCGA and 'gtex' for GTEx data
    :return: sample of the kind requested
    :rtype: str
    """
    for s in samples:
        if kind == 'tumor':
            if s.startswith('TCGA') and s.endswith('01'):
                yield s
        elif kind == 'normal':
            if s.startswith('TCGA') and s.endswith('11'):
                yield s
        elif kind == 'gtex':
            if s.startswith('GTEX'):
                yield s
        else:
            raise RuntimeError('Invalid argument for "kind". Valid options are "tumor", "normal", "gtex". '
                               'User submitted: {}'.format(kind))


def load_gene_map():
    """
    Dictionary mapping gene ID to gene name

    :return: Gene map
    :rtype: dict
    """
    return pickle.load(open(os.path.join(_cwd, 'gene_map.pickle'), 'rb'))


def load_mab_targets():
    """
    Returns sorted list of MAB cancer drug targets

    :return: Sorted gene list
    :rtype: list
    """
    path = os.path.join(_cwd, 'cancer-MAB-gene-targets.txt')
    return sorted([x.strip() for x in open(path, 'r').readlines()])


def load_drug_gene_tissue_table():
    """
    Loads Dataframe containing Drug, Gene, Tissue, and DE information

    :return: Dataframe
    :rtype: pd.DataFrame
    """
    return pd.read_csv(os.path.join(_cwd, 'drugs', 'drug-gene-tissue-deseq2.tsv'), sep='\t', index_col=0)


def load_gene_tissue_drug_map():
    """
    Tissue and Drug information for a given gene

    :return: Mapping
    :rtype: dict(str, list(str))
    """
    df = load_drug_gene_tissue_table()
    return {k: {'tissue': set(g['tissue'].tolist()), 'drug': set(g['generic_name'].tolist())}
            for k, g in df.groupby('gene')}


def load_ucsf_genes():
    """
    Returns sorted list of UCSF genes

    :return: Sorted gene list
    :rtype: list
    """
    path = os.path.join(_cwd, 'UCSF-genes.csv')
    return sorted([x.strip() for x in open(path, 'r').readlines()])


def load_civic_genes():
    """
    Returns sorted list of genes from CIViC

    :return: Sorted gene list
    :rtype: list
    """
    path = os.path.join(_cwd, 'data/civic-genes.txt')
    return sorted([x.strip() for x in open(path, 'r').readlines()])


def load_cosmic_dataframe():
    """
    Returns cosmic dataframe

    :return: Dataframe of COSMIC genes
    :rtype: pd.DataFrame
    """
    return pd.read_csv(os.path.join(_cwd, 'data/census_all_1-26-2018.tsv'))
