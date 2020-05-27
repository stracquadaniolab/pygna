import logging
import networkx as nx
import numpy as np
import scipy
import sys
import matplotlib.pyplot as plt
import pygna.reading_class as rc
import pygna.diagnostic as diag
import multiprocessing
import time


class KSTest:
    def __init__(self, distribution, network, diz={}):
        self.__distribution = distribution

        self.__network = network
        self.__diz = diz

        # print(type(self.__network))
        if (type(self.__network) is nx.Graph) or (type(self.__network) is nx.DiGraph):
            self.__universe = set(self.__network.nodes())
        elif type(self.__network) is dict:
            self.__universe = set(self.__network.keys())
            # print (self.__network [list(self.__network.keys())[0]])
        else:
            logging.error("Unknown network type: %s" % type(self.__network))
            sys.exit(-1)

    def apply_test(self, geneset):
        # mapping geneset
        mapped_geneset = sorted(list(set(geneset).intersection(self.__universe)))
        others = list(self.__universe.difference(set(mapped_geneset)))
        if mapped_geneset == 0:
            return 0, 0
        else:
            logging.info(
                "Mapped %d genes out of %d." % (len(mapped_geneset), len(geneset))
            )
            observed = self.__distribution(
                self.__network, mapped_geneset, self.__diz, observed_flag=True
            )
            baseline = self.__distribution(
                self.__network, others, self.__diz, observed_flag=True
            )
            logging.info("Observed Distribution evaluated")
            # iterations, the null distribution is a vector where ones represent rejected hypothesis
            stats, pvalue = scipy.stats.ks_2samp(observed, baseline, alternative='greater')
            # computing empirical pvalue

            return stats, pvalue, len(mapped_geneset), len(geneset)


###############################################################################
###  TEST STATISTICS FOR SINGLE GENESET  ######################################
###############################################################################


def degree_distribution(network, geneset, diz={}, observed_flag=False):
    """ degree distribution of the geneset """
    degree = nx.degree(network)
    distribution = np.array([degree[g] for g in geneset])
    return distribution




################################################################################
######### Degree Distribution test #############################################
################################################################################

# FIXME: problem with scipy version 1.2.1 , does not have alternative for test
def test_degree_distribution(
    network_file: "network file",
    geneset_file: "GMT geneset file",
    output_table: "output results table, use .csv extension",
    setname: "Geneset to analyse" = None,
    size_cut: "removes all genesets with a mapped length < size_cut" = 20,
    results_figure: "barplot of results, use pdf or png extension" = None,
    ):
    """
        Performs degree distribution test.
        Kolmogorov-Smirnov statistic on 2 samples.
        H0 is that the geneset is drawn from the same distribution of all the other nodes.
        H0 rejected if statistic is greater.
    """
    network = rc.ReadTsv(network_file)

    geneset = rc.ReadGmt(geneset_file).get_geneset(setname)

    setnames = [key for key in geneset.keys()]

    output1 = out.Output(
        network_file, output_table, "test_degree", geneset_file, setnames
    )
    logging.info("Results file = " + output1.output_table_results)
    output1.create_st_table_empirical()

    st_test = KS.KSTest(KS.degree_distribution, network)

    for setname, item in geneset.items():

        item = set(item)
        if len(item) > size_cut:
            observed, pvalue, n_mapped, n_geneset = st_test.apply_test(item)
            logging.info("Setname:" + setname)
            if n_mapped < size_cut:
                logging.info(
                    "%s remove from results since nodes mapped are < %d"
                    % (setname, size_cut)
                )
            else:
                logging.info("Observed: %g p-value: %g" % (observed, pvalue))

                output1.update_st_table_empirical(
                    setname, n_mapped, n_geneset, 1, observed, pvalue, 0, 0
                )

    output1.close_temporary_table()
    if results_figure:
        paint.paint_datasets_stats(output1.output_table_results, results_figure, alternative='greater')



#TODO: painter class to be removed
class Painter:
    def __init__(self, network, output_folder, name, diz={}):

        self.network = network
        self.output_folder = output_folder
        self.name = name
        self.diz = diz

        if (type(self.network) is nx.Graph) or (type(self.network) is nx.DiGraph):
            self.universe = set(self.network.nodes())
        elif type(self.network) is dict:
            self.universe = set(self.network.keys())
        else:
            logging.error("Unknown network type: %s" % type(self.network))
            sys.exit(-1)


class Painter_RW(Painter):
    def __init__(self, network, output_folder, name, RW_dict, geneset, diz={}):
        Painter.__init__(self, network, output_folder, name, diz={})
        self.RW_nodes = list(RW_dict["nodes"])
        self.geneset = list(geneset)
        self.geneset_index = [list(RW_dict["nodes"]).index(i) for i in geneset]
        self.RW_matrix = RW_dict["matrix"][self.geneset_index, :][:, self.geneset_index]

    def plot_matrix(self, show_labels=False):

        logging.info("Plotting figure as " + str(self.output_folder + self.name))

        fig, axes = plt.subplots(1, 1)
        fig.subplots_adjust(left=0.2, right=0.99, bottom=0.2, top=0.99)
        axes.set_title("Diffusion matrix of geneset's nodes")
        mask = np.eye(self.RW_matrix.shape[0])

        logging.info(
            str(np.array(self.geneset)[:, np.newaxis].shape) + str(self.RW_matrix.shape)
        )

        if show_labels:
            tab_conversion = pd.read_table(
                "/home/viola/Desktop/geneset-network-analysis/primary_data/entrez_name.tsv",
                sep="\t",
            )
            labels = []
            for i in self.geneset:
                name = tab_conversion[tab_conversion["EntrezID"] == int(i)][
                    "Symbol"
                ].values.tolist()
                if len(name) > 0:
                    labels.append(str(name[0]))
                else:
                    labels.append(i)

            g2 = sns.heatmap(
                self.RW_matrix - mask,
                cmap="OrRd",
                square=True,
                ax=axes,
                xticklabels=labels,
                yticklabels=labels,
                cbar=True,
                vmin=0,
                vmax=np.max(self.RW_matrix - mask),
            )
            g2.set_yticklabels(g2.get_yticklabels(), rotation=0, fontsize=6)
            g2.set_xticklabels(g2.get_xticklabels(), rotation=90, fontsize=6)
        else:
            g2 = sns.heatmap(self.RW_matrix, cmap="OrRd", mask=mask, ax=axes, cbar=True)

        fig.savefig(
            self.output_folder + "RW_matrix_" + self.name + ".pdf", format="pdf"
        )

    def restrict_network(self, geneset):

        self.__network = nx.DiGraph(
            self.__network.subgraph(max(nx.connected_components(self.__network), key=len))
        )

    def add_annotation_nodes_in_geneset(self, geneset, annotation="in_subset"):

        logging.info("adding annotation to nodes in geneset")
        logging.info("%d genes in geneset,%d mapped in network ",
            (len(geneset), len(self.__universe.intersection(set(geneset)))),
        )
        geneset = self.__universe.intersection(set(geneset))
        dict_nodes = {}
        for n in self.__network.nodes():
            if n in geneset[setname]:
                dict_nodes[n] = True
            else:
                dict_nodes[n] = False
        nx.set_node_attributes(self.__network, annotation, dict_nodes)

    def add_weights_edges(self, dict_edges):
        pass

    def draw_graphml(self):

        """
        Draws a .graphml file from the network
        """
        # mapping geneset
        logging.info("Drawing graph for " + str(name))
        nx.write_graphml(
            self.__network, self.__output_folder + self.__name + ".graphml"
        )
