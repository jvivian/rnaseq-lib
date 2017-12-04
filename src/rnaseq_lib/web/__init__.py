import urllib2

import mygene
import requests
import xmltodict
from bs4 import BeautifulSoup
from progressbar import ProgressBar

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
        return None
    else:
        return r
