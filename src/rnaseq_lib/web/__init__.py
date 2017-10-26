import time
import urllib2
import pandas as pd

import mygene
import requests
import xmltodict
from bs4 import BeautifulSoup
from collections import defaultdict
from progressbar import ProgressBar

from rnaseq_lib.utils import rexpando

bar = ProgressBar()


def get_drug_target_from_wiki(drug):
    """
    Scrape wikipedia for the target of a drug

    :param str drug: Drug to lookup
    :return: Drug target
    :rtype: str
    """
    # Look for wiki page
    url = 'https://en.wikipedia.org/wiki/'
    try:
        page = urllib2.urlopen(url + drug)
    except urllib2.HTTPError:
        print 'Page not found for: {}'.format(drug)
        return None

    # Parse page
    soup = BeautifulSoup(page, 'html.parser')

    # Look for table via HTML tags
    name_box = soup.find('table', {'class': 'infobox'})
    if name_box:

        # Look for "Target", next item in the list should be the drug name
        name = name_box.text.strip()
        if 'Target' in name:
            return name.split('\n')[name.split('\n').index('Target') + 1]
        else:
            print '{} has no listed Target'.format(drug)
            return None

    else:
        print 'No table found for {}'.format(drug)
        return None


def get_info_from_wiki(drug):
    """
    Scrape wikipedia for all information about a drug

    :param str drug: Drug to lookup
    :return: Information on the wikipedia page
    :rtype: str
    """
    # Look for wiki page
    url = 'https://en.wikipedia.org/wiki/'
    try:
        page = urllib2.urlopen(url + drug)
    except urllib2.HTTPError:
        print 'Page not found for: {}'.format(drug)
        return None

    # Parse page
    soup = BeautifulSoup(page, 'html.parser')

    # Scrape all content
    name_box = soup.find('div', {'class': 'mw-content-ltr'})
    if name_box:
        return name_box.text.strip()
    else:
        print 'No table found for {}'.format(drug)
        return None


def openfda_query_drug(drug):
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


def openfda_get_drugs_by_query(query, limit=100, field='indications_and_usage'):
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


def openfda_drugs_to_dataframe(drugs):
    """
    Convert a list of drugs to an OpenFDA DataFrame

    :param list(str) drugs:
    :return: DataFrame of OpenFDA information
    """
    # Create dictionary to store table info
    info = defaultdict(list)

    # For each drug, check openFDA for info
    bar = ProgressBar()
    for drug in bar(drugs):
        drug = drug.lower()
        r = openfda_query_drug(drug)
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


def get_drug_usage_nih(drug):
    """
    Gets drug uasage information from NIH API

    :param str drug:
    :return: Usage section from NIH
    :rtype: str
    """
    # Make request
    params = {'drug_name': drug}
    url = 'https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json'
    r = _rget(url, params=params)
    if r:

        # Get "set ID" to query SPLS for detailed info
        setid = None
        for data in r.json()['data']:
            if drug in data['title'] or drug.upper() in data['title']:
                print 'Found set ID for: {}'.format(data['title'])
                setid = data['setid']
        if setid:

            # Make request
            url = 'https://dailymed.nlm.nih.gov/dailymed/services/v2/spls/{}.xml'.format(setid)
            r = _rget(url)
            if r:

                # I hate XML with all my being
                xml = xmltodict.parse(r.content)
                comp = xml['document']['component']['structuredBody']['component']

                # Look for usage tag
                content = None
                for sec in comp:
                    if 'USAGE' in sec['section']['code']['@displayName']:
                        content = str(sec['section'])

                # Parse
                if content:
                    remove = ['[', '(', ')', ']', "u'", "'", '#text', 'OrderedDict',
                              'u@styleCode', 'uitalics', 'ulinkHtml', 'href', ',']
                    for item in remove:
                        content = content.replace(item, '')
                    return content.replace('displayName', '\n\n')

                else:
                    print 'No content found'
                    return None

            else:
                return None

        else:
            print 'No set ID found for {}'.format(drug)
            return None

    else:
        return None


def find_gene_given_alias(alias, strict=True):
    """
    Queries MyGene to look for valid gene aliases

    :param str alias: Gene alias/name to query
    :param bool strict: If strict, will only return gene if part of gene_map set
    :return: gene
    :rtype: str
    """
    from rnaseq_lib.tissues import get_gene_map

    # Create valid gene set from gene_map
    gene_map = get_gene_map()
    valid_genes = set(gene_map.keys() + gene_map.values())

    # MyGene query
    mg = mygene.MyGeneInfo()
    try:
        hits = mg.query(alias, fields='symbol,ensembl.gene')['hits']
    except KeyError:
        hits = []

    # Iterate through hits for gene
    gene = None
    for hit in hits:
        if hit['symbol'] in valid_genes or hit['symbol'].upper() in valid_genes:
            gene = hit['symbol']
            break

        # If no matching symbol is found, look for ensemble name
        else:
            try:
                if hit['ensembl']['gene'] in valid_genes:
                    gene = hit['ensembl']['gene']
            except KeyError:
                pass

    if gene in valid_genes:
        return gene
    elif gene and not strict:
        print 'WARNING: Gene found, but invalid: {}'.format(gene)
    return gene


def _rget(url, params=None):
    """
    Wrapper for requests.get that checks status code

    :param str url: Request URL
    :param dict params: Parameters for request
    :return: Request from URL or None if status code != 200
    :rtype: requests.models.Response
    """
    r = requests.get(url, params=params)
    if r.status_code != 200:
        # print 'Error Status Code {}'.format(r.status_code)
        return None
    else:
        return r
