import pandas as pd
import progressbar
import requests

bar = progressbar.ProgressBar()


def create_civic_drug_disease_dataframe():
    """
    Hits CIViC's API to build a dataframe on genes, cancers, and drugs

    :return: Dataframe of CIVIC information
    :rtype: pd.DataFrame
    """
    payload = {'count': 999}  # Only 295 genes at time of writing
    request = requests.get('https://civic.genome.wustl.edu/api/genes/', params=payload)
    assert request.status_code == 200, 'Request failed: {}'.format(request.status_code)

    records = request.json()['records']

    # Create records in the format:  Cancer, Gene, Drugs, Variant-Name, Aliases, Description
    info = []
    for record in bar(records):
        gene = record['name']
        description = record['description']
        aliases = ','.join(record['aliases'])
        variants = record['variants']

        for variant in variants:
            var_name = variant['name']
            var_id = variant['id']
            var = requests.get('https://civic.genome.wustl.edu/api/variants/{}'.format(var_id)).json()
            for evidence in var['evidence_items']:
                disease = evidence['disease']['name']
                drugs = ','.join([x['name'] for x in evidence['drugs']])
                info.append((disease, drugs, gene, aliases, var_name, description))

    # Create DataFrame
    labels = ['Cancer', 'Drugs', 'Gene', 'Aliases', 'Variant-Name', 'Description']
    return pd.DataFrame.from_records(info, columns=labels)
