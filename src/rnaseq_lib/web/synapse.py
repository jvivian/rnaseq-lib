import os

from synapseclient import Synapse, File

expression = 'syn11311347'
metadata = 'syn11311931'
objects = 'syn11515014'
pathway = 'syn11585314'


def upload_file(file_path, login, parent, description=None):
    """
    Uploads file to Synapse. Password must be stored in environment variable SYNAPSE_PASS

    :param str file_path: Path to file
    :param str login: Login (usually an email address)
    :param str parent: Parent Synapse ID (example: syn12312415) where file will be placed
    :param str description: Optional description to add
    """
    description = '' if None else description
    f = File(file_path, description=description, parent=parent)

    syn = _syn_login(login)
    syn.store(f)


def download_file(syn_id, login, download_location='.'):
    """
    Synapse ID of file to download. Password must be stored in environment variable SYNAPSE_PASS

    :param str syn_id: Synapse ID
    :param str login: Login (usually an email address)
    :param str download_location: Download location for file
    """
    syn = _syn_login(login)
    syn.get(syn_id, downloadLocation=download_location)


def _syn_login(login):
    """
    Login to synapse. Password must be stored in environment variable SYNAPSE_PASS

    :param str login:
    :return: Synapse instance
    :rtype: instance
    """
    assert 'SYNAPSE_PASS' in os.environ, 'SYNAPSE_PASS must be set as an environment variable'
    syn = Synapse()
    syn.login(login, os.environ['SYNAPSE_PASS'])
    return syn


# TODO: Finish eventually
def mirror_synapse_dir(syn_id, login):
    """
    Mirror a synapse directory's contents. Password must be stored in environment variable SYNAPSE_PASS

    :param str syn_id: Synapse ID
    :param str login: Login (usually an email address)
    """
    pass
