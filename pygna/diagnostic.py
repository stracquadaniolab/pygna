import matplotlib.pyplot as plt
import seaborn as sns
import logging
import numpy as np
import networkx as nx


def plot_degree(degree_object, output_file):

    """
    Diagnosis tool for the degree object
    """

    D = dict(degree_object)
    degrees = {k: v for k, v in D.items()}
    degree_values=np.array(list(degrees.values()))

    fig, axes = plt.subplots(1, figsize=(10, 10))
    g1 = sns.distplot(degree_values, hist=True, ax=axes)
    
    key_max = max(degrees.keys(), key=(lambda k: degrees[k]))
    g1 = sns.distplot([degrees[key_max]], hist=False, kde=False, rug=True, color='r', ax=axes)
    axes.annotate('%s: %d' %(key_max, degrees[key_max]), xy=(degrees[key_max], 0),
                xytext=(degrees[key_max], axes.dataLim.y1/2),
                arrowprops=dict(arrowstyle="->")
                )    

    g1 = sns.distplot([np.median(degree_values)], hist=False, kde=False, rug=True, color='r', ax=axes)
    axes.annotate('median %f' %np.median(degree_values), xy=(np.median(degree_values), 0),
                            xytext=(np.median(degree_values), axes.dataLim.y1/2),
                            arrowprops=dict(arrowstyle="->")
                            )    

    sns.despine(ax=axes, top=True, bottom=False, right=True, left=True)
    g1.set_ylabel("Density")
    g1.set_xlabel("Node Degree")

    if output_file.endswith('.pdf'):
        plt.savefig(output_file, format="pdf")
    elif output_file.endswith('.png'):
        plt.savefig(output_file, format="png")
    else:
        logging.warning('The null distribution figure can only be saved in pdf or png, forced to png')
        fig.savefig(output_file+'.png', format="png")

def plot_connected_components(c_components, output_file):

    """
    Diagnosis tool for the connected components object.

    Creates the histogram of the components length, to analyse the relationship between the lcc
    and the other c_components, and prints some overall stats about the connected components
    """
    c_components_len = [len(k) for k in c_components]

    logging.info("Number of cc %d" % (len(c_components_len)))
    logging.info("First five cc" + str(c_components_len[0:5]))
    logging.info("Mean length of cc %d" % (np.mean(c_components_len)))

    fig, axes = plt.subplots(1, figsize=(10, 10))
    g1 = sns.distplot(c_components_len, hist=True, kde=False,ax=axes,norm_hist=False)
    g1 = sns.distplot([np.max(c_components_len)], hist=False, kde=False, rug=True, color='r', ax=axes,norm_hist=False)

    axes.annotate('LCC: %d' %np.max(c_components_len), xy=(np.max(c_components_len), 0),
                xytext=(np.max(c_components_len)-10,axes.dataLim.y1/4), arrowprops=dict(arrowstyle="->"))    

    sns.despine(ax=axes, top=True, bottom=False, right=True, left=True)
    g1.set_ylabel("Number of CC")
    g1.set_xlabel("Size of CC")

    if output_file.endswith('.pdf'):
        plt.savefig(output_file, format="pdf")
    elif output_file.endswith('.png'):
        plt.savefig(output_file, format="png")
    else:
        logging.warning('The null distribution figure can only be saved in pdf or png, forced to png')
        fig.savefig(output_file+'.png', format="png")

def plot_diffusion_matrix(nodes, matrix, filename, show_labels=False):

    """
    Diagnosis tool for a diffusion matrix.

    Shows the weighted adjacency matrix that is the output of a build process
    """

    logging.info("Plotting figure as " + str(filename))
    fig, axes = plt.subplots(1)
    axes.imshow(matrix, cmap="PuBu")
    if show_labels == True:
        pass
    plt.show()
    fig.savefig(filename + ".pdf", format="pdf")


def plot_null_distribution(null_distribution, observed, output_file, setname):

    """
    Saves the density plot of the null distribution and pinpoints the observed value
    """

    fig, axes = plt.subplots(1, figsize=(10, 10))
    g1 = sns.distplot(null_distribution, hist=True, kde=True, rug=False, ax=axes)
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

    logging.info(
        "Output for diagnostic null distribution: " + output_file
    )
    if output_file.endswith('.pdf'):
        fig.savefig(output_file, format="pdf")
    elif output_file.endswith('.png'):
        fig.savefig(output_file, format="png")
    else:
        logging.warning('The null distribution figure can only be saved in pdf or png, forced to png')
        fig.savefig(output_file+'.png', format="png")


