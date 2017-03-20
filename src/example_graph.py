#!/usr/bin/python

# Load modules.
import math
import sys, argparse

from hhio import save_index_gene, save_edge_list, save_gene_score

# Parse arguments.
def get_parser():
    description = 'Generate an example vertex-weighted graph.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-igf', '--index_gene_file', type=str, required=True, help='Index-gene filename')
    parser.add_argument('-elf', '--edge_list_file', type=str, required=True, help='Edge list filename')
    parser.add_argument('-gsfa', '--gene_score_file_a', type=str, required=True, help='Gene-score filename')
    parser.add_argument('-gsfb', '--gene_score_file_b', type=str, required=True, help='Gene-score filename')
    parser.add_argument('-pf', '--plot_file', type=str, required=False, help='Plot filename')
    return parser

# Run script.
def run(args):
    # Initialize nodes, etc.
    node_list = list()
    edge_list = list()
    score_list_a = list()
    score_list_b = list()
    pos = dict()

    # Add a clique.
    letters = 'abcdefgh'
    node_list += [letters[i] for i in range(8)]
    for i in range(8):
        for j in range(i+1, 8):
            edge_list += [(letters[i], letters[j])]
    score_list_a += [5]*8
    score_list_b += [0.1]*8
    for i in range(8):
        pos[letters[i]] = (math.cos(2*math.pi*i/8)-2, math.sin(2*math.pi*i/8))


    # Add an approximate clique.
    letters = 'ijklmnopij'
    node_list += [letters[i] for i in range(8)]
    for i in range(8):
        edge_list += [(letters[i], letters[i+1])]
        edge_list += [(letters[i], letters[i+2])]
    score_list_a += [10]*8
    score_list_b += [0.4, 0.4, 0.5, 0.6, 0.6, 0.6, 0.5, 0.4]
    for i in range(8):
        pos[letters[i]] = (math.cos(2*math.pi*i/8)+2, math.sin(2*math.pi*i/8))

    # Add a cycle.
    letters = 'qrstuvwxq'
    node_list += [letters[i] for i in range(8)]
    for i in range(8):
        edge_list += [(letters[i], letters[i+1])]
    score_list_a += [3]*8
    score_list_b += [0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2]
    for i in range(8):
        pos[letters[i]] = (math.cos(2*math.pi*i/8), math.sin(2*math.pi*i/8)-2)

    # Add a point connecting the clique, approximate clique, and cycle, add an isolated point.
    node_list += ['y']
    edge_list += [('a', 'y'), ('m', 'y'), ('s', 'y')]
    score_list_a += [5]
    score_list_b += [1.0]
    pos['y'] = (0.0, 0.0)

    # Save Hierarchical HotNet input.
    index_to_gene = dict((i+1, gene) for i, gene in enumerate(node_list))
    gene_to_index = dict((gene, i) for i, gene in index_to_gene.items())
    gene_to_score_a = dict((gene, score) for gene, score in zip(node_list, score_list_a))
    gene_to_score_b = dict((gene, score) for gene, score in zip(node_list, score_list_b))

    save_index_gene(args.index_gene_file, index_to_gene)
    save_edge_list(args.edge_list_file, edge_list, gene_to_index)
    save_gene_score(args.gene_score_file_a, gene_to_score_a)
    save_gene_score(args.gene_score_file_b, gene_to_score_b)

    # Draw the graph for illustration.
    if args.plot_file:
        import networkx as nx
        import matplotlib.pyplot as plt

        G = nx.Graph()
        G.add_nodes_from(node_list)
        G.add_edges_from(edge_list)

        plt.figure(figsize=(15, 5))
        plt.subplot(121)
        node_labels = dict((v, r'${}$'.format(v)) for v in node_list)
        nx.draw_networkx_nodes(G, pos=pos, nodelist=node_list, node_color=score_list_a, alpha=0.2, cmap=plt.cm.YlOrRd)
        nx.draw_networkx_edges(G, pos=pos, alpha=0.2)
        nx.draw_networkx_labels(G, pos=pos, labels=node_labels, alpha=0.0)
        plt.axis('off')

        plt.subplot(122)
        node_labels = dict((v, r'${}$'.format(v)) for v in node_list)
        nx.draw_networkx_nodes(G, pos=pos, nodelist=node_list, node_color=score_list_b, alpha=0.2, cmap=plt.cm.YlOrRd)
        nx.draw_networkx_edges(G, pos=pos, alpha=0.2)
        nx.draw_networkx_labels(G, pos=pos, labels=node_labels, alpha=0.0)
        plt.axis('off')

        plt.savefig(args.plot_file, bbox_inches='tight')

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
