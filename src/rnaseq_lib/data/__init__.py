import os

import pandas as pd

_cwd = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


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


def get_ucsf_genes():
    """
    Returns sorted list of UCSF genes

    :return: Sorted gene list
    :rtype: list
    """
    path = os.path.join(_cwd, 'UCSF-genes.csv')
    return sorted([x.strip() for x in open(path, 'r').readlines()])


def get_civic_genes():
    """
    Returns sorted list of genes from CIViC

    :return: Sorted gene list
    :rtype: list
    """
    path = os.path.join(_cwd, 'data/civic-genes.txt')
    return sorted([x.strip() for x in open(path, 'r').readlines()])
