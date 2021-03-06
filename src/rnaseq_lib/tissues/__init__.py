import os
import re
from collections import defaultdict

import pandas as pd
from rnaseq_lib.data import load_gene_map
from rnaseq_lib.gtf import get_protein_coding_genes
from rnaseq_lib.utils import flatten

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

# dtype for metadata-expression dataframe
dtype = {'id': object,
         'reads': int,
         'size_MB': float,
         'sex': object,
         'tissue': object,
         'seq_site': object,
         'weight': float,
         'height': float,
         'mapped_reads': float,
         'race': object,
         'age': float,
         'qc': float,
         'dataset': object,
         'tumor': object,
         'type': object,
         'label': object}


def identify_tissue_from_str(content):
    """
    Identifies possible tissue(s) referenced by a given string by matching terms that are "complete words"

    :param str content: Text to examine for terms associated with tissues
    :return: Possible tissues referenced in input string
    :rtype: set(str)
    """
    td_map = tissue_disease_mapping()

    # \b Matches the empty string, but only at the beginning or end of a word
    return set([k for k, v in td_map.iteritems() if
                any([term for term in v if re.search(r'\b{}\b'.format(term), content, flags=re.IGNORECASE)])])


def subset_by_dataset(df):
    """
    Returns tumor, normal, and GTEx dataframe subsets

    :param pd.Dataframe df: Dataframe containing a "dataset" and "tumor" column
    :return: One dataframe per dataset
    :rtype: tuple(pd.Dataframe, pd.Dataframe, pd.Dataframe)
    """
    tumor = df[df.label == 'tcga-tumor']
    normal = df[df.label == 'tcga-normal']
    gtex = df[(df.label == 'gtex')]
    return tumor, normal, gtex


def tissue_disease_mapping():
    """
    Maps tissue types to words associated with cancers of that tissue

    :return: Tissue / disease term mapping
    :rtype: dict(str, list(str))
    """
    return {
        'Adrenal': ['adrenal', 'adrenocortical', 'cortical', 'oncocytic', 'myxoid'],
        'Bladder': ['bladder', 'urothelial'],
        'Blood': ['blood', 'leukemia', 'myeloma', 'hematologic'],
        'Bone': ['bone', 'osteosarcoma', 'ewing'],
        'Brain': ['brain', 'astrocytoma', 'neurocytoma', 'choroid', 'plexus', 'neuroepithelial',
                  'ependymal', 'Pleomorphic Xanthoastrocytoma', 'fibrillary', 'giant-cell', 'glioblastoma',
                  'multiforme', 'gliomatosis', 'cerebri', 'gliosarcoma',
                  'hemangiopericytoma', 'medulloblastoma', 'medulloepithelioma', 'meningeal', 'neuroblastoma',
                  'neurocytoma', 'oligoastrocytoma', 'optic', 'ependymoma', 'pilocytic', 'pinealoblastoma',
                  'pineocytoma', 'meningioma', 'subependymoma', 'retinoblastoma', 'neuro', 'glioma', 'ganglioglioma',
                  'oligodendroglioma', 'ganglioglioma'],
        'Breast': ['breast'],
        'Cervix': ['cervix', 'cervical'],
        'Colon': ['colon', 'rectal', 'colorectal', 'intestine', 'intestinal', 'bowel', 'lynch', 'rectum'],
        'Esophagus': ['esophagus', 'barrett', 'esophageal', 'oropharyngeal', 'oropharynx', 'laryngeal'],
        'Kidney': ['kidney', 'renal ', 'nephron', 'nephroma', 'wilm', 'chromophobe', 'rhabdoid', 'papillary renal'],
        'Liver': ['liver', 'hepatic', 'Hepatocellular', 'parenchymal', 'cholang'],
        'Lung': ['lung', 'small-cell carcinoma', 'non-small-cell carcinoma', 'small cell carcinoma',
                 'non small cell carcinoma', 'non small-cell carcinoma', 'mesothelioma'],
        'Ovary': ['ovary', 'ovarian', 'endometrioid', 'fallopian', 'cord', 'bronchogenic'],
        'Pancreas': ['pancreas', 'pancreatic', 'cystadenocarcinomas'],
        'Prostate': ['prostate'],
        'Skin-Head': ['head', 'neck', 'skin', 'basal', 'melanoma', 'oral', 'merkel', 'salivary'],
        'Stomach': ['stomach', 'gastric', 'bile', 'cholangiocarcinoma', 'Gastrointestinal'],
        'Testis': ['testis', 'testicular', 'testes', 'gonad', ],
        'Thyroid': ['thyroid', 'papillary carcinoma'],
        'Uterus': ['uterus', 'uterine', 'endometrial', 'ureteral', 'gestational']
    }


def grep_cancer_terms(content, replace_newlines_with_periods=True, comprehensive=False):
    """
    Returns sentences with cancer terms

    :param str content: String containing sentences to check for cancer terms
    :param bool replace_newlines_with_periods: If True, replaces newlines with periods so they count as "sentences"
    :param bool comprehensive: if True, adds all values from tissue_disease_mapping
    :return: Sentences with matches
    :rtype: list(str)
    """
    terms = {'cancer', 'leukemia', 'carcinoma', 'squamous', 'lymphoma', 'malignant',
             'metastasis', 'metastatic', 'sarcoma', 'tumor', 'adenocarcinoma'}

    # Add all terms from tissue_disease_mapping to grep list
    terms = terms.union(set(flatten(tissue_disease_mapping().values()))) if comprehensive else terms

    # Replace newlines with periods
    content = content.replace('\n', '.') if replace_newlines_with_periods else content

    # Return sentences that include terms
    return [x for x in content.split('.') if any(y for y in terms if y.upper() in x.upper())]


def validate_genes(input_genes):
    from rnaseq_lib.web import find_gene_given_alias

    # Check if user submits single gene as a string
    input_genes = [input_genes] if type(input_genes) == str else input_genes

    # Create valid_genes list
    gene_map = load_gene_map()
    valid_genes = set(gene_map.keys() + gene_map.values())

    # Iterate through genes
    genes = []
    for g in input_genes:

        # Check if listed Gene is valid
        gene = None
        if g in valid_genes:
            gene = g

        # If gene invalid, search MyGene
        if not gene:
            gene = find_gene_given_alias(g)

        # Append gene or None
        genes.append(gene)

    return genes


def subset_protein_coding_genes(df, gtf_path):
    """
    Return dataframe with protein-coding genes

    :param pd.DataFrame df: Genes by samples gene expression dataframe
    :param str gtf_path:
    :return: Subset dataframe
    :rtype: pd.DataFrame
    """
    genes = get_protein_coding_genes(gtf_path)
    return df.loc[genes]


def parse_pathway(filepath):
    """
    Parses source -> targets TSV and returns mapping dictionaries

    :param str filepath: File path to pathway TSV
    :return:  Source and Target mappings
    :rtype: tuple(dict(str, list), dict(str, list))
    """
    source_targets = {}
    target_sources = defaultdict(list)

    # Parse pathway file to generate first dict
    with open(filepath, 'r') as f:
        f.readline()
        for line in f:
            line = line.split()
            tf, targets = line[0], line[1:]
            source_targets[tf] = targets

    for tf, targets in source_targets.iteritems():
        for target in targets:
            target_sources[target].append(tf)

    return source_targets, target_sources
