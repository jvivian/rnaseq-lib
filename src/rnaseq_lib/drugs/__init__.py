import os

import pandas as pd

__location__ = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__))))


def get_anova():
    return pd.read_excel(os.path.join(__location__, 'data/cancerrxgene/ANOVA_OUTPUT.xlsx'))


def get_cell_line():
    return pd.read_excel(os.path.join(__location__, 'data/cancerrxgene/Cell_Lines_Details.xlsx'))


def get_drugs():
    return pd.read_excel(os.path.join(__location__, 'data/cancerrxgene/Screened_Compounds.xlsx'))


def get_dose_response():
    return pd.read_excel(os.path.join(__location__, 'data/cancerrxgene/Screened_Compounds.xlsx'))
