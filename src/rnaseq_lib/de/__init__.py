import numpy as np

def get_l2fc(a, b):
    return np.log2(b + 1) - np.log2(a + 1)
