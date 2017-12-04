import os

import pandas as pd

from rnaseq_lib.web import _rget


# Gene Functions
def find_gene(query, database='genes', human_only=True):
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
    label = find_gene(gene)
    r = _get(label)
    return _parse_gene_get(r.text)


# Drug Functions
def find_drug_label(drug):
    r = _kegg_search(operation='find', database='drug', query=drug)
    for line in r.text.split('\n'):
        try:
            return line.split()[0]
        except IndexError:
            print 'No drug found, check spelling'


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
        print 'Response found, but no drug class provided'
    else:
        return values


# Internal Functions
def _kegg_search(operation, database, query):
    url = 'http://rest.kegg.jp'
    return _rget(os.path.join(url, operation, database, query))


def _get(query):
    return _kegg_search(operation='get', database='', query=query)


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
