import os

from synapseclient import Synapse, File

expression = 'syn11311347'


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

    assert 'SYNAPSE_PASS' in os.environ, 'SYNAPSE_PASS must be set as an environment variable'

    syn = Synapse()
    syn.login(login, os.environ)
    syn.store(f)
