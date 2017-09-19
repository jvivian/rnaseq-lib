import urllib2

import requests
import xmltodict
from bs4 import BeautifulSoup


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

    # Look for table
    name_box = soup.find('table', {'class': 'infobox'})

    if not name_box:
        print 'No table found for {}'.format(drug)

    # Look for Target, next item should be the Drug
    name = name_box.text.strip()
    if 'Target' in name:
        return name.split('\n')[name.split('\n').index('Target') + 1]
    else:
        print '{} has no listed Target'.format(drug)
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


def grep_cancer_terms(content):
    """
    Returns sentences with cancer terms

    :param str content: String containing sentences to check for cancer terms
    :return: Sentences with matches
    :rtype: list(str)
    """
    terms = {'cancer', 'leukemia', 'carcinoma', 'squamous', 'lymphoma',
             'malignant', 'metastasis', 'metastatic', 'sarcoma', 'tumor'}

    return [x.upper() for x in content.split('.') if any(y.upper() for y in terms if y in x)]


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
        print 'Error Status Code {}'.format(r.status_code)
        return None
    else:
        return r
