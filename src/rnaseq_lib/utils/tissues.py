import os
import pickle

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


def return_samples():
    """
    Returns sample dictionary which maps TCGA and GTEx samples to a tissue.
    Synapse ID: syn10296681

    :return: Tissues are keys are list of samples are values
    :rtype: dict[list]
    """
    return pickle.load(open(os.path.join(__location__, 'samples.pickle'), 'rb'))
