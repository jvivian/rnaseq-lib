import os
import pickle
import pandas as pd

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


def return_samples():
    """
    Returns sample dictionary which maps TCGA and GTEx samples to a tissue.
    Synapse ID: syn10296681

    :return: Tissues are keys are list of samples are values
    :rtype: dict[list]
    """
    return pickle.load(open(os.path.join(__location__, 'data/samples.pickle'), 'rb'))


def get_gene_map():
    """
    Dictionary mapping gene ID to gene name

    :return: Gene map
    :rtype: dict
    """
    return pickle.load(open(os.path.join(__location__, 'data/gene_map.pickle'), 'rb'))


def map_genes(df):
    gene_map = get_gene_map()
    genes = [gene_map[x.split('.')[0]] for x in df.index]
    df.index = genes
    return df


def get_mab_targets():
    """
    Returns sorted list of MAB cancer drug targets

    :return: Sorted gene list
    :rtype: list
    """
    path = os.path.join(__location__, 'data/cancer-MAB-gene-targets.txt')
    return sorted([x.strip() for x in open(path, 'r').readlines()])


def get_ucsf_genes():
    """
    Returns sorted list of UCSF genes

    :return: Sorted gene list
    :rtype: list
    """
    path = os.path.join(__location__, 'data/UCSF-genes.csv')
    return sorted([x.strip() for x in open(path, 'r').readlines()])


def get_ucsf_subset(df):
    """
    Subset UCSF dataframe and return.

    :param pd.DataFrame df: Input Dataframe in the format of "Genes by Samples"
    :return: Subset of Dataframe that only includes UCSF genes
    :rtype: pd.DataFrame
    """
    df = map_genes(df)

    ucsf_genes = get_ucsf_genes()
    ucsf_genes = [x for x in ucsf_genes if x in df.index]

    return df.loc[ucsf_genes]
