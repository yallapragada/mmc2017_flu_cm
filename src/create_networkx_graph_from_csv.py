import networkx as nx
import pandas as pd
import sys
import os


graph = nx.Graph()


# create pandas dataframe from csv
def create_pd_df_from_csv(filename):
    csv_df = pd.read_csv(filename)
    return csv_df


def check_if_node_exists(node):
    if node in graph:
        return True
    else:
        return False


# create nodes for each unique source & target residues in csv
def create_nodes(csv_df):

    source_nodes = csv_df.x.unique()
    target_nodes = csv_df.y.unique()

    for node in source_nodes:
        node_str = str(node).upper()
        tokens = node_str.split('_')
        protein = tokens[0]
        residue_number = tokens[1]

        # don't create a node if one exists already
        if check_if_node_exists(node_str) is False:
            graph.add_node(node_str, protein=protein)
        else:
            print('node exists for %s %s' % (protein, residue_number))

    for node in target_nodes:
        node_str = str(node).upper()
        tokens = node_str.split('_')
        protein = tokens[0]
        residue_number = tokens[1]

        # don't create a node if one exists already
        if check_if_node_exists(node_str) is False:
            graph.add_node(node_str, protein=protein)
        else:
            print('node exists for %s %s' % (protein, residue_number))


def create_edges(csv_df):

    for index, row in csv_df.iterrows():
        source = str(row['x']).upper()
        target = str(row['y']).upper()
        weight = row['weight']
        graph.add_edge(source, target, weight=weight)


def load(filename):
    csv_df = create_pd_df_from_csv(filename)
    create_nodes(csv_df=csv_df)
    create_edges(csv_df=csv_df)


if __name__ == '__main__':
    file = sys.argv[1]
    load(file)
    out_file = file.split('.')[0] + '.graphml'
    nx.write_graphml(graph, out_file)