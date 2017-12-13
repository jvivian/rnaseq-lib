import io
import os

import pandas as pd
from Bio.KEGG.KGML import KGML_parser

from rnaseq_lib.web import _rget


# Gene Functions
def find_gene_label(query, database='genes', human_only=True):
    r = _kegg_search(operation='find', database=database, query=query)
    if r:
        for hit in r.text.split('\n'):
            hit = hit.split('\t')[0]
            if human_only:
                if hit.startswith('hsa'):
                    return str(hit)
            else:
                return str(hit)
    else:
        return None


def get_gene_dataframe(gene):
    label = find_gene_label(gene)
    r = _get(label)
    return _parse_gene_get(r.text)


def get_gene_names_from_label(label):
    genes = set()
    r = _get(label)
    for line in r.text.split('\n'):
        if line.startswith('NAME'):
            line = line.split()[1:]
            genes.add(line[0].rstrip(','))
    return genes


# Drug Functions
def find_drug_label(drug):
    r = _kegg_search(operation='find', database='drug', query=drug)
    for line in r.text.split('\n'):
        try:
            return line.split()[0]
        except IndexError:
            print 'No entry for {}, check spelling'.format(drug)


def get_drug_info(drug):
    label = find_drug_label(drug)
    r = _get(label)
    return r.text


def get_drug_class(drug):
    text = get_drug_info(drug)
    values = []

    # Flag for once the BRITE is reached
    brite = False
    for line in text.split('\n'):
        line = line.strip().split()
        if not line:
            break
        if line[0] == 'BRITE':
            brite = True
        if line[0].startswith('L') and brite:
            values.append(' '.join(line))

    # Handle returning values (or notify user if no values found)
    if text and not values:
        print 'Response found for {}, but no drug class provided'.format(drug)
    else:
        return values


def get_drug_disease(drug):
    text = get_drug_info(drug)
    values = []

    disease = True
    for line in text.split('\n'):
        line = line.strip().split()
        if not line:
            break
        if line[0].isupper() and disease:
            break
        if line[0].upper() == 'DISEASE':
            disease = True
        if disease:
            line.remove('DISEASE')
            values.append(' '.join(line))


# Pathway Functions
def find_pathway(query):
    return _kegg_search(operation='find', database='pathway', query=query)


def get_kgml(pathway):
    kgml = _get(pathway, form='kgml').text

    # Wrap text in a file handle for KGML parser
    f = io.BytesIO(kgml.encode('utf-8'))
    k = KGML_parser.read(f)
    return k


def get_pathway_info(pathway):
    return _get(pathway).text


def pathway_image(pathway):
    from IPython.display import Image
    k = get_kgml(pathway)
    return Image(k.image)


def get_genes_from_pathway(pathway):
    kgml = _get(pathway, form='kgml').text

    # Wrap text in a file handle for KGML parser
    f = io.BytesIO(kgml.encode('utf-8'))
    k = KGML_parser.read(f)

    genes = set()
    for gene in k.genes:
        g = get_gene_names_from_label(gene)
        genes = genes.union(g)

    return genes


# Internal Functions
def _kegg_search(operation, database, query=None, form=None):
    # Set arguments to empty strings if None
    query = '' if query is None else query
    form = '' if form is None else form

    # Define base URL
    url = 'http://rest.kegg.jp'

    # Make get request
    request = os.path.join(url, operation, database, query, form)
    return _rget(request)


def _get(query, form=None):
    return _kegg_search(operation='get', database='', query=query, form=form)


def _parse_drug_line(line, index=1):
    drug_value = None
    for item in line[index:]:
        if item.startswith('D') and not item.startswith('DG'):
            drug_value = item
            break
    return drug_value


def _parse_gene_get(text):
    # Entries in the text we want to include
    entries = ['ORTHOLOGY', 'ORGANISM', 'PATHWAY', 'NETWORK', 'DISEASE']

    # Create first entry (header) in dataframe
    records = [('Category', 'Kegg', 'Entry')]

    # Flag variable for specifying current entry
    current_entry = None

    for line in text.split('\n'):

        # If line is empty or contains trailing slashes from request
        if line.isspace() or line.startswith('/'):
            break

        # Split line by whitespace as there is no standard delimiter
        line = line.split()

        # First line
        if line[0] == 'ENTRY':
            pass

        # Name row has no Kegg value
        elif line[0] == 'NAME':
            records.append(('Name', None, ' '.join(line[1:])))

        # Definition also has no Kegg value
        elif line[0] == 'DEFINITION':
            records.append(('Definition', None, ' '.join(line[1:])))

        # Drug target row is non-standard, requires separate parsing
        elif line[0] == 'DRUG_TARGET':
            drug_value = _parse_drug_line(line, index=2)
            records.append(('Drug', drug_value, line[1]))
            current_entry = 'Drug'

        # Continuation of drug entries
        elif current_entry == 'Drug':
            drug_value = _parse_drug_line(line, index=1)
            records.append(('Drug', drug_value, line[0]))

        # First line for standard entries
        elif line[0] in entries:
            name = line[0].lower().capitalize()
            records.append((name, line[1], ' '.join(line[1:])))
            current_entry = name

        # Stopping condition -- the rest of the table is hard to parse and less interesting
        elif line[0] == 'BRITE':
            break

        # If in an entry and not yet at the stopping condition
        elif current_entry:
            records.append((current_entry, line[0], ' '.join(line[1:])))

    # Return pandas dataframe
    return pd.DataFrame.from_records(records)
