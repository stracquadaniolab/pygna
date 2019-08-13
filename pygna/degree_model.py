import networkx as nx
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import logging
import pandas as pd
import seaborn as sns
from pygna import output
from pygna.utils import YamlConfig


class DegreeModel(object):
    def __init__(self, network_prob=0.5, vip_prob=1, n_nodes=10, vip_percentage=0.1):

        self.n_nodes = n_nodes
        self.graph = nx.Graph()
        self.network_prob = network_prob
        self.vip_prob = vip_prob
        self.n_vip = math.ceil(vip_percentage * self.n_nodes)
        self.nodes = ["N" + str(i) for i in range(n_nodes)]
        self.cluster_dict = {}

    def set_nodes(self, nodes_names):
        self.nodes = nodes_names
        self.n_nodes = len(nodes_names)

    def create_graph(self):

        reject = True
        logging.info("Reject=" + str(reject))
        while reject:
            graph = generate_graph_vip(
                self.n_nodes,
                self.n_vip,
                network_prob=self.network_prob,
                vip_prob=self.vip_prob,
                node_names=self.nodes,
            )
            LCC = max(nx.connected_components(graph), key=len)
            reject = len(LCC) != self.n_nodes
            logging.info("Reject=" + str(reject))
            logging.info("Nodes: %d, in LCC: %d" % (self.n_nodes, len(LCC)))

        self.graph = graph

    def plot_graph(self):
        pass

    def write_network(self, output_file):

        self.network_file = output_file

        logging.info("Network written on %s" % (output_file))

        if output_file.endswith(".tsv"):
            nx.write_edgelist(self.graph, output_file, data=False, delimiter="\t")
        else:
            logging.error("output file format unknown")

    def write_genelist(self, output_file):

        self.genelist_file = output_file

        clusters = nx.get_node_attributes(self.graph, "cluster")

        for i in set(clusters.values()):
            c = "cluster_" + str(i)
            self.cluster_dict[c] = {}
            self.cluster_dict[c]["descriptor"] = "cluster"
            self.cluster_dict[c]["genes"] = [
                str(j) for j in clusters.keys() if clusters[j] == i
            ]

        if output_file.endswith(".gmt"):
            output.print_GMT(self.cluster_dict, self.genelist_file)
        else:
            logging.error("output file format unknown")


def generate_graph_vip(n_nodes, n_vip, network_prob=0.5, vip_prob=1, node_names=None):
    """
    This function creates a graph with n_nodes number of vertices and a matrix
    block_model that describes the intra e inter- block connectivity.
    The nodes_in_block is parameter, list, to control the number of nodes in each cluster
    """

    if not node_names:
        node_names = range(n_nodes)

    edges = []
    G = nx.Graph()

    list_temp = [(n_nodes - n_vip) * [0]]
    list_temp.append(n_vip * [1])

    cluster = np.array([val for sublist in list_temp for val in sublist])
    np.random.shuffle(cluster)

    prob = [network_prob, vip_prob]
    p = [prob[i] for i in cluster]

    for i in range(n_nodes):
        G.add_node(node_names[i], cluster=cluster[i], prob=p[i])

    for i in range(n_nodes):
        for j in range(i + 1, n_nodes):
            if np.random.binomial(1, p[i]) == 1:
                edges.append((node_names[i], node_names[j]))

    G.add_edges_from(edges)
    return G


def plot_vip_graph(graph, output_folder=None):

    nodes = graph.nodes()
    colors = ["#b15928", "#1f78b4"]
    cluster = nx.get_node_attributes(graph, "cluster")
    labels = [colors[cluster[n]] for n in nodes]
    layout = nx.spring_layout(graph)

    plt.figure(figsize=(13.5, 5))
    plt.subplot(1, 3, 1)
    nx.draw(
        graph,
        nodelist=nodes,
        pos=layout,
        node_color="#636363",
        node_size=50,
        edge_color="#bdbdbd",
    )
    plt.title("Observed network")

    plt.subplot(1, 3, 2)
    plt.imshow(nx.adjacency_matrix(graph).toarray(), cmap="OrRd")
    plt.title("Adjacency Matrix")

    plt.subplot(1, 3, 3)
    legend = []
    for ix, c in enumerate(colors):
        legend.append(mpatches.Patch(color=c, label="C%d" % ix))

    nx.draw(
        graph,
        nodelist=nodes,
        pos=layout,
        node_color=labels,
        node_size=50,
        edge_color="#bdbdbd",
    )
    plt.legend(handles=legend, ncol=len(colors), mode="expand", borderaxespad=0)
    plt.title("VIP clustering")

    plt.savefig(output_folder + "VIP.pdf", bbox_inches="tight")


def plot_adjacency(graph, output_folder, prefix):

    plt.figure(figsize=(13.5, 5))

    plt.subplot(1, 1, 1)
    plt.imshow(nx.adjacency_matrix(graph).toarray(), cmap="OrRd")
    plt.title("Adjacency Matrix")

    plt.savefig(output_folder + prefix + "VIP.png")


def generate_hdn_network(
    output_folder,
    prefix,
    n_nodes: "number of nodes in the network" = 1000,
    network_prob: "probability of connection in the network" = 0.005,
    vip_prob: "probability of the connection of VIP" = 0.3,
    vip_percentage=0.05,
    number_of_simulations=5,
    ):

    """ This function generates a simulated network using the VIP model
    """

    dm = DegreeModel(
        network_prob=network_prob,
        vip_prob=vip_prob,
        n_nodes=n_nodes,
        vip_percentage=vip_percentage,
    )

    for i in range(number_of_simulations):
        dm.create_graph()
        dm.write_network(output_folder + prefix + "_s_" + str(i) + "_network.tsv")
        dm.write_genelist(output_folder + prefix + "_s_" + str(i) + "_genes.gmt")
        plot_adjacency(dm.graph, output_folder, prefix=prefix + "_s_" + str(i))
