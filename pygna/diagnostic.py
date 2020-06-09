import matplotlib.pyplot as plt
import seaborn as sns
import logging
import numpy as np
from matplotlib.offsetbox import AnchoredText
import networkx as nx


def plot_degree(degree_object: nx.Graph, output_file: str):
    """
    Diagnosis tool for the degree object

    :param degree_object: the graph to plot
    :param output_file: the path to save the file
    """

    D = dict(degree_object)
    degrees = {k: v for k, v in D.items()}
    degree_values = np.array(list(degrees.values()))

    fig, axes = plt.subplots(1, figsize=(10, 10))
    g1 = sns.distplot(degree_values, hist=True, ax=axes)

    key_max = max(degrees.keys(), key=(lambda k: degrees[k]))
    g1 = sns.distplot([degrees[key_max]], hist=False, kde=False, rug=True, color='r', ax=axes)
    axes.annotate('%s: %d' % (key_max, degrees[key_max]), xy=(degrees[key_max], 0),
                  xytext=(degrees[key_max], axes.dataLim.y1 / 2), arrowprops=dict(arrowstyle="->"))

    g1 = sns.distplot([np.median(degree_values)], hist=False, kde=False, rug=True, color='r', ax=axes)
    axes.annotate('median %f' % np.median(degree_values), xy=(np.median(degree_values), 0),
                  xytext=(np.median(degree_values), axes.dataLim.y1 / 2), arrowprops=dict(arrowstyle="->"))

    sns.despine(ax=axes, top=True, bottom=False, right=True, left=True)
    g1.set_ylabel("Density")
    g1.set_xlabel("Node Degree")

    if output_file.endswith('.pdf'):
        plt.savefig(output_file, format="pdf")
    elif output_file.endswith('.png'):
        plt.savefig(output_file, format="png")
    else:
        logging.warning('The null distribution figure can only be saved in pdf or png, forced to png')
        fig.savefig(output_file + '.png', format="png")


def plot_connected_components(c_components: nx.connected_components, output_file: str) -> None:
    """
    Diagnosis tool for the connected components object.
    Creates the histogram of the components length, to analyse the relationship between the lcc
    and the other c_components, and prints some overall stats about the connected components

    :param c_components: the list of the connected components
    :param output_file:  the path to save the file
    """
    c_components_len = [len(k) for k in c_components]

    logging.info("Number of cc %d" % (len(c_components_len)))
    logging.info("First five cc" + str(c_components_len[0:5]))
    logging.info("Mean length of cc %d" % (np.mean(c_components_len)))

    fig, axes = plt.subplots(1, figsize=(10, 10))
    g1 = sns.distplot(c_components_len, hist=True, kde=False, ax=axes, norm_hist=False)
    g1 = sns.distplot([np.max(c_components_len)], hist=False, kde=False, rug=True, color='r', ax=axes, norm_hist=False)

    axes.annotate('LCC: %d' % np.max(c_components_len), xy=(np.max(c_components_len), 0),
                  xytext=(np.max(c_components_len) - 10, axes.dataLim.y1 / 4), arrowprops=dict(arrowstyle="->"))

    sns.despine(ax=axes, top=True, bottom=False, right=True, left=True)
    g1.set_ylabel("Number of CC")
    g1.set_xlabel("Size of CC")

    if output_file.endswith('.pdf'):
        plt.savefig(output_file, format="pdf")
    elif output_file.endswith('.png'):
        plt.savefig(output_file, format="png")
    else:
        logging.warning('The null distribution figure can only be saved in pdf or png, forced to png')
        fig.savefig(output_file + '.png', format="png")


def plot_diffusion_matrix(nodes: list, matrix: np.matrix, filename: str, show_labels: bool = False) -> None:
    """
    Diagnosis tool for a diffusion matrix.
    Shows the weighted adjacency matrix that is the output of a build process

    :param nodes: the network nodes
    :param matrix: the diffusion matrix
    :param filename: the path to save the file
    :param show_labels: if labels should be plotted
    """

    logging.info("Plotting figure as " + str(filename))
    fig, axes = plt.subplots(1)
    axes.imshow(matrix, cmap="PuBu")
    if show_labels:
        pass
    plt.show()
    fig.savefig(filename + ".pdf", format="pdf")


def plot_null_distribution(null_distribution: list, observed: list, output_file: str, setname: str,
                           alternative: str = "greater") -> None:
    """
    Saves the density plot of the null distribution and pinpoints the observed value

    :param null_distribution: the list with the values from the null distribution
    :param observed: list of the observed genes
    :param output_file: the path to save the file
    :param setname: the name of the gene set
    :param alternative: use "greater" if you want to take the genes with greater than the observed value
    """

    fig, axes = plt.subplots(1, figsize=(8, 6))
    g1 = sns.distplot(null_distribution, hist=True, kde=True, rug=False, ax=axes)
    if alternative == "greater":
        if len(null_distribution[null_distribution > observed]):
            g3 = sns.distplot(null_distribution[null_distribution > observed], hist=False, kde=False, rug=True,
                              rug_kws={'height': 1 / 50}, color="r", ax=axes)
    else:
        if len(null_distribution[null_distribution < observed]):
            g3 = sns.distplot(null_distribution[null_distribution < observed], hist=False, kde=False, rug=True,
                              rug_kws={'height': 1 / 50}, color="r", ax=axes)
    ymax = axes.dataLim.y1
    xmax = axes.dataLim.x1
    print('xmax %f' % xmax)
    g4 = axes.stem([observed], [ymax / 2], "r", "r--")

    sns.despine(ax=axes, top=True, bottom=False, right=True, left=True)
    anchored_text = AnchoredText("Observed:%1.1E" % observed, loc=1, prop={'fontsize': 12, 'color': 'r'},
                                 **{'frameon': False})
    axes.add_artist(anchored_text)
    axes.set_xlabel('Statistics', fontsize=12)
    axes.set_ylabel('Density', fontsize=12)
    logging.info("Output for diagnostic null distribution: " + output_file)
    if output_file.endswith('.pdf'):
        fig.savefig(output_file, format="pdf")
    elif output_file.endswith('.png'):
        fig.savefig(output_file, format="png")
    else:
        logging.warning('The null distribution figure can only be saved in pdf or png, forced to png')
        fig.savefig(output_file + '.png', format="png")
