import matplotlib.pyplot as plt
import seaborn as sns
import logging
import numpy as np
import networkx as nx


def connected_components_diagnostic(c_components):
    """
    Diagnosis tool for the connected components object.

    Creates the histogram of the components length, to analyse the relationship between the lcc
    and the other c_components, and prints some overall stats about th connected components
    """
    c_components_len = [len(k) for k in c_components]

    logging.info("Number of cc %d" % (len(c_components)))
    logging.info("First five cc" + str(c_components_len[0:5]))
    logging.info("Mean length of cc %d" % (np.mean(c_components_len)))

    color = "#f00000"

    plt.figure()
    sns.distplot(c_components_len, kde=False, color="b")
    plt.show()

    # print(cc_len[0:5])
    # print("median cc %g " % (np.min(cc_len)))
    # print("Min cc: %g Max cc: %g " % (np.min(cc_len), np.max(cc_len)))


def plot_diffusion_matrix(nodes, matrix, filename, show_labels=False):

    """
    Diagnosis tool for a diffusion matrix.

    Shows the weighted adjacency matrix that is the output of a build process
    """

    print(matrix)
    logging.info("Plotting figure as " + str(filename))
    fig, axes = plt.subplots(1)
    axes.imshow(matrix, cmap="PuBu")
    if show_labels == True:
        pass
    plt.show()
    fig.savefig(filename + ".pdf", format="pdf")

    # print(cc_len[0:5])
    # print("median cc %g " % (np.min(cc_len)))
    # print("Min cc: %g Max cc: %g " % (np.min(cc_len), np.max(cc_len)))


def plot_null_distribution(null_distribution, observed, output_folder, setname=None):

    fig, axes = plt.subplots(1, figsize=(10, 10))
    g1 = sns.distplot(null_distribution, hist=True, kde=True, rug=False, ax=axes)
    # g2= sns.distplot(null_distribution, kde=False, rug=False, color= "r", ax=axes)
    if len(null_distribution[null_distribution > observed]):
        g3 = sns.distplot(
            null_distribution[null_distribution > observed],
            hist=False,
            kde=False,
            rug=True,
            color="r",
            ax=axes,
        )
    g4 = axes.stem([observed], [0.01], "r", "r--")

    sns.despine(ax=axes, top=True, bottom=False, right=True, left=True)
    axes.annotate(
        "observed=%d" % observed,
        xy=(observed, 0.01),
        xytext=(observed - 45, 0.01),
        color="r",
        fontsize=10,
    )
    # plt.show()

    logging.info(
        "output file for figure: " + output_folder + setname + "_null_distribution.pdf"
    )
    fig.savefig(output_folder + setname + "_null_distribution.pdf", format="pdf")


def draw_graph(
    network: "network file", output_file: "name graphml network for visualisation"
):
    """
        Plotting graph
    """

    for setname in geneset:
        print(setname)
        new_network = nx.Graph()
        new_network_minimal = nx.Graph()
        for s in geneset[setname]:
            if s in network.nodes():
                l0 = np.inf
                for t in geneset[setname]:
                    if (t in network.nodes()) & (s != t):
                        path = nx.shortest_path(network, source=s, target=t)
                        if l0 > len(path):
                            l0 = len(path)
                            new_network_minimal.add_path(path)
                        new_network.add_path(path)

        dict_nodes = {}
        for n in new_network.nodes():
            if n in geneset[setname]:
                dict_nodes[n] = True
            else:
                dict_nodes[n] = False
        nx.set_node_attributes(new_network, "in_subset", dict_nodes)
        nx.write_graphml(new_network, output_folder + setname + ".graphml")
