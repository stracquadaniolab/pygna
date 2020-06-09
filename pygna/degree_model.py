import networkx as nx
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import logging
from pygna import output


class DegreeModel(object):
    def __init__(self, network_prob: float = 0.5, vip_prob: float = 1, n_nodes: int = 10,
                 vip_percentage: float = 0.1):
        self.n_nodes = n_nodes
        self.graph = nx.Graph()
        self.network_prob = network_prob
        self.vip_prob = vip_prob
        self.n_vip = math.ceil(vip_percentage * self.n_nodes)
        self.nodes = ["N" + str(i) for i in range(n_nodes)]
        self.cluster_dict = {}

    def set_nodes(self, nodes_names: list) -> None:
        """
        Set the name of the nodes

        :param nodes_names: the list with the nodes name
        """
        self.nodes = nodes_names
        self.n_nodes = len(nodes_names)

    def create_graph(self) -> None:
        """
        Create a graph from the nodes
        """
        reject = True
        logging.info("Reject=" + str(reject))
        while reject:
            graph = generate_graph_vip(self.n_nodes, self.n_vip, network_prob=self.network_prob, vip_prob=self.vip_prob,
                                       node_names=self.nodes)
            LCC = max(nx.connected_components(graph), key=len)
            reject = len(LCC) != self.n_nodes
            logging.info("Reject=" + str(reject))
            logging.info("Nodes: %d, in LCC: %d" % (self.n_nodes, len(LCC)))

        self.graph = graph

    def plot_graph(self):
        # Todo
        pass

    def write_network(self, output_file: str) -> None:
        """
        Write on file the network as an edge list

        :param output_file: the file path where to save the network
        """
        self.network_file = output_file

        logging.info("Network written on %s" % output_file)

        if output_file.endswith(".tsv"):
            nx.write_edgelist(self.graph, output_file, data=False, delimiter="\t")
        else:
            logging.error("output file format unknown")

    def write_genelist(self, output_file: str) -> None:
        """
        Write the GMT gene list on file

        :param output_file: the file path where to save the gene list
        """
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


def generate_graph_vip(n_nodes: int, n_vip: int, network_prob: float = 0.5, vip_prob: float = 1,
                       node_names: list = None) -> nx.Graph:
    """
    This function creates a graph with n_nodes number of vertices and a matrix block_model that describes the intra e inter-block connectivity.
    The nodes_in_block is parameter, list, to control the number of nodes in each cluster

    :param n_nodes: number of nodes in the network
    :param n_vip: number of VIP to create
    :param network_prob: probability of connection in the network
    :param vip_prob: probability of connection of the vip
    :param node_names: list of nodes for the network
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


def plot_vip_graph(graph: nx.Graph, output_folder: str = None) -> None:
    """
    Plot the VIP graph on the specific folder

    :param graph: the graph to plot
    :param output_folder: the folder path where to save the file
    """
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


def plot_adjacency(graph: nx.Graph, output_folder: str, prefix: str) -> None:
    """
    Plot the adjacency matrix on file

    :param graph: the graph to plot
    :param output_folder: the folder where to save the file
    :param prefix: the prefix to give to the file
    """
    plt.figure(figsize=(13.5, 5))

    plt.subplot(1, 1, 1)
    plt.imshow(nx.adjacency_matrix(graph).toarray(), cmap="OrRd")
    plt.title("Adjacency Matrix")

    plt.savefig(output_folder + prefix + "VIP.png")


def generate_hdn_network(output_folder: str, prefix: str, n_nodes: int = 1000, network_prob: float = 0.005,
                         hdn_probability: float = 0.3, hdn_percentage: float = 0.05, number_of_simulations: int = 5):
    """
    This function generates a simulated network using the VIP model

    :param output_folder: the output folder path
    :param prefix: the prefix of the file to be saved
    :param n_nodes: the number of nodes in the network
    :param network_prob: probability of connection in the network
    :param hdn_probability: probability of connection of the VIP
    :param hdn_percentage: percentage of connection
    :param number_of_simulations: number of simulation to be performed
    """

    dm = DegreeModel(
        network_prob=network_prob,
        vip_prob=hdn_probability,
        n_nodes=n_nodes,
        vip_percentage=hdn_percentage,
    )

    for i in range(number_of_simulations):
        dm.create_graph()
        dm.write_network(output_folder + prefix + "_s_" + str(i) + "_network.tsv")
        dm.write_genelist(output_folder + prefix + "_s_" + str(i) + "_genes.gmt")
        plot_adjacency(dm.graph, output_folder, prefix=prefix + "_s_" + str(i))
