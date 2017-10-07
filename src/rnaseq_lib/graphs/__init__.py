import os
from itertools import combinations

import pandas as pd


def output_gephi(df, node_attrs, output_dir):
    """
    Output two CSVs used by Gephi to create graphs:
        nodes.csv
        edges.csv

    :param pd.DataFrame df: Dataframe of input data
    :param list(str) node_attrs: List of attributes in DF to create nodes from
    :param str output_dir: Path to output directory
    """
    # Create nodes
    nodes = []
    for attr in node_attrs:
        nodes.extend(map(lambda x: (x, attr), df[attr].tolist()))
    nodes = pd.DataFrame.from_records(nodes, columns=('Label', 'Type'))

    # Create edges
    edges = set()
    for row in df.iterrows():
        i, row = row
        for j, k in combinations(node_attrs, 2):
            edges.add((int(nodes[nodes.Label == row[j]].index[0]), int(nodes[nodes.Label == row[k]].index[0])))
    edges = pd.DataFrame.from_records(list(edges), columns=('Source', 'Target'))

    # Output
    nodes.to_csv(os.path.join(output_dir, 'nodes.csv'), index=True, index_label='Id')
    edges.to_csv(os.path.join(output_dir, 'edges.csv'), index=False)
