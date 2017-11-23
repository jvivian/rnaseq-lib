import time
from collections import defaultdict

import pandas as pd
from progressbar import ProgressBar
from rnaseq_lib.utils import rexpando
from rnaseq_lib.web import _rget
from rnaseq_lib.tissues import grep_cancer_terms


def query_drug(drug):
    """
    Search OpenFDA drug label API by name

    :param str drug: Drug to query by
    :return: Response containing first match
    :rtype: requests.models.Response
    """
    # Clean input
    drug = drug.lower()

    # Look for drug label on openFDA
    for d in [drug, drug.capitalize(), drug.upper()]:
        for name in ['generic_name', 'brand_name']:
            url = 'https://api.fda.gov/drug/label.json?search=openfda.{}:{}&limit=1'.format(name, d)
            time.sleep(0.1)
            r = _rget(url)
            if r:
                return r
    return None


def drugs_by_query(query, limit=100, field='indications_and_usage'):
    """
    Search OpenFDA API for drugs that match query

    :param str query: Query string to search by
    :param int limit: Limit request, cannot be higher than 100
    :param str field: OpenFDA field to search. Example: "openfda.brand_name"
    :return: Return drugs by query term
    :rtype: list(tuple(str, str))
    """
    assert limit <= 100, 'OpenFDA API does not allow a limit higher than 100'
    url = 'https://api.fda.gov/drug/label.json?search={}:{}&limit={}'.format(field, query, limit)
    r = _rget(url)
    if r:
        # Convert to JSON
        r = r.json()

        # Get total number of hits
        total_terms = r['meta']['results']['total']
        print 'Found a total of {} terms'.format(total_terms)

        # Collect first batch
        drugs = [(','.join(x['openfda']['brand_name']), ','.join(x['openfda']['generic_name']))
                 for x in r['results']]

        # Collect results in batches
        bar = ProgressBar()
        skip = limit
        for _ in bar(xrange(int(total_terms) / limit)):
            time.sleep(1)
            print 'Collecting samples {} - {}'.format(skip, skip + limit)
            r = _rget(url + '&skip={}'.format(skip))
            if r:
                r = r.json()
                drugs.extend([(','.join(x['openfda']['brand_name']), ','.join(x['openfda']['generic_name']))
                              for x in r['results']])
                skip += limit
        return drugs
    else:
        return None


def drugs_to_dataframe(drugs):
    """
    Convert a list of drugs to an OpenFDA DataFrame

    :param list(str) drugs:
    :return: DataFrame of OpenFDA information
    """
    # Create dictionary to store table info
    info = defaultdict(list)

    # For each drug, query openFDA
    bar = ProgressBar()
    for drug in bar(drugs):
        drug = drug.lower()
        r = query_drug(drug)
        if r:
            hits = rexpando(r.json()).results

            # If more than one hit is returned, find exact match
            if len(hits) != 1:
                for h in hits:
                    if drug in h.openfda.brand_name.lower() or drug in h.openfda.generic_name.lower():
                        res = h
            else:
                res = hits[0]

            # Collect info if description references cancer
            if 'indications_and_usage' in res:
                usage = res.indications_and_usage[0] if type(
                    res.indications_and_usage) is list else res.indications_and_usage
                if [x for x in usage.replace('\n', '.').split('.') if 'cancer' in x]:

                    # Get generic name
                    try:
                        info['generic_name'].append(str(res.openfda.generic_name).strip("u'[]"))
                    except AttributeError:
                        info['generic_name'].append(None)

                    # Get brand name
                    try:
                        info['brand_name'].append(str(res.openfda.brand_name).strip("u'[]"))
                    except AttributeError:
                        info['brand_name'].append(None)

                    # Get usage
                    try:
                        info['usage'].append(str(res.indications_and_usage).strip('[]'))
                    except AttributeError:
                        info['usage'].append(None)

                    # Get mechanism of action
                    try:
                        info['mech_action'].append(str(res.mechanism_of_action).strip('[]'))
                    except AttributeError:
                        try:
                            info['mech_action'].append(str(res.clinical_pharmacology).strip('[]'))
                        except AttributeError:
                            info['mech_action'].append(None)

    return pd.DataFrame.from_dict(info)


def indications_and_usage(drug):
    """
    Returns indications and usage, split by cancer terms

    :param str drug: Name of drug
    :return: Indications and usage split by cancer terms
    :rtype: str
    """
    r = query_drug(drug)
    if not r:
        print 'No drug found'
        return None

    return '\n\n'.join(grep_cancer_terms(str(rexpando(
        r.json()['results'][0])['indications_and_usage'])))