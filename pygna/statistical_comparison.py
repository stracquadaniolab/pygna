import logging
import networkx as nx
import numpy as np
import sys
import multiprocessing
import pygna.statistical_test as st


class StatisticalComparison:
    """
    This class implements the statistical analysis comparison between two genesets.
    Please refer to the single method documentation for the returning values
    """

    def __init__(self, comparison_statistic, network: nx.Graph, n_proc: int = 1, diz: dict = {}):
        """
        :param comparison_statistic: the method do be applied to perform the statistical analysis
        :param network: the network to be used during the analysis
        :param n_proc: the number of CPU used to perform the elaboration
        :param diz: the dictionary used for the genesets
        """
        self.__comparison_statistic = comparison_statistic
        self.__network = network
        self.__diz = diz
        self.__n_proc = n_proc

        if (type(self.__network) is nx.Graph) or (type(self.__network) is nx.DiGraph):
            self.__universe = set(self.__network.nodes())
        elif type(self.__network) is dict:
            self.__universe = set(self.__network.keys())

        else:
            logging.error("Unknown network type: %s" % type(self.__network))
            sys.exit(-1)

    def comparison_empirical_pvalue(self, genesetA: set, genesetB: set, alternative: str = "less", max_iter: int = 100,
                                    keep: bool = False) -> [int, float, float, int, int]:
        """
        Calculate the empirical value between two genesets

        :param genesetA: the first geneset to compare
        :param genesetB: the second geneset to compare
        :param alternative: the pvalue selection of the observed genes
        :param max_iter: the maximum number of iterations
        :param keep: if the geneset B should not be kept
        :return observed, pvalue, null_distribution, len(mapped_genesetA), len(mapped_genesetB): the list with the data calculated
        """

        # mapping genesets
        mapped_genesetA = sorted(list(set(genesetA) & self.__universe))
        mapped_genesetB = sorted(list(set(genesetB) & self.__universe))
        logging.info("Mapped %d genes out of %d from genesetA" % (len(mapped_genesetA), len(genesetA)))
        logging.info("Mapped %d genes out of %d from genesetB" % (len(mapped_genesetB), len(genesetB)))

        # Observed statistical values
        observed = self.__comparison_statistic(self.__network, mapped_genesetA, mapped_genesetB, self.__diz)

        # iterations
        null_distribution = StatisticalComparison.get_comparison_null_distribution_mp(self, mapped_genesetA,
                                                                                      mapped_genesetB, max_iter, keep)
        pvalue = 1
        if alternative == "greater":
            pvalue = (np.sum(null_distribution >= observed) + 1) / (float(len(null_distribution) + 1))
        else:
            pvalue = (np.sum(null_distribution <= observed) + 1) / (float(len(null_distribution) + 1))

        return observed, pvalue, null_distribution, len(mapped_genesetA), len(mapped_genesetB)

    def get_comparison_null_distribution_mp(self, genesetA: list, genesetB: list, max_iter: int = 100,
                                            keep: bool = False) -> np.ndarray:
        """
        Calculate the null distribution between two genesets with multiple CPUs

        :param genesetA: the first geneset to compare
        :param genesetB: the second geneset to compare
        :param max_iter: maximum number of iteration to perform
        :param keep: if the geneset B should not be kept
        :return: the array with null distribution
        """

        n_trial = int(max_iter / self.__n_proc)
        logging.info("n_proc = %d, each computing %d permutations " % (int(self.__n_proc), n_trial))

        # Don't set up the multicore architecture if only one core is given
        if self.__n_proc == 1:
            null_distribution = StatisticalComparison.get_comparison_null_distribution(self, genesetA, genesetB,
                                                                                       max_iter, keep)
        else:
            p = multiprocessing.Pool(self.__n_proc)
            results = [p.apply_async(StatisticalComparison.get_comparison_null_distribution,
                                     args=(self, genesetA, genesetB, n_trial, keep)) for w in
                       list(range(1, self.__n_proc + 1))]

            null_distribution = np.array([])
            for r in results:
                null_distribution = np.hstack((null_distribution, np.array(r.get())))
            p.close()

        return np.asarray(null_distribution)

    def get_comparison_null_distribution(self, genesetA: list, genesetB: list, n_samples: int, keep: bool) -> list:
        """
        Calculate the null distribution between two genesets with single CPU

        :param genesetA: the first geneset to compare
        :param genesetB: the second geneset to compare
        :param n_samples: the number of samples to be taken
        :param keep: if the geneset B should not be kept
        :return: the random distribution calculated
        """
        np.random.seed()
        random_dist = []

        # association statistic, B kept, A bootstrapped
        if keep:
            for i in range(n_samples):
                random_sample_A = np.random.choice(list(self.__universe), len(genesetA), replace=False)
                random_dist.append(self.__comparison_statistic(self.__network, set(random_sample_A), set(genesetB),
                                                               self.__diz))
        # association statistic, B kept, A bootstrapped
        else:
            for i in range(n_samples):
                random_sample_A = np.random.choice(list(self.__universe), len(genesetA), replace=False)
                random_sample_B = np.random.choice(list(self.__universe), len(genesetB), replace=False)
                random_dist.append(self.__comparison_statistic(self.__network, set(random_sample_A),
                                                               set(random_sample_B), self.__diz))
        return random_dist


###############################################################################
###  TEST STATISTICS FOR COMPARISONS  #########################################
###############################################################################

def comparison_shortest_path(network: nx.Graph, genesetA: list, genesetB: list, diz: dict) -> float:
    """
    Evaluate the shortest path between two genesets

    :param network: the graph representing the network
    :param genesetA: the first geneset list
    :param genesetB: the second geneset list
    :param diz: the dictionary containing the nodes name and index
    """
    n = np.array([diz["nodes"].index(i) for i in genesetA])
    m = np.array([diz["nodes"].index(i) for i in genesetB])
    if len(n) == 0 or len(m) == 0:
        logging.info("Geneset length is equal to 0")
        sys.exit()
    cum_sum = calculate_sum(n, m, diz)
    cum_sum += calculate_sum(m, n, diz)
    d_AB = cum_sum / float(len(genesetA) + len(genesetB))
    d_A = st.geneset_localisation_statistic(network, genesetA, diz)
    d_B = st.geneset_localisation_statistic(network, genesetB, diz)
    return d_AB - (d_A + d_B) / 2


def calculate_sum(n: np.ndarray, m: np.ndarray, diz: dict) -> np.ndarray:
    """
    Evaluate the sum of the columns of two matrices

    :param n: the first column
    :param m: the second column
    :param diz: the dictionary containing the data
    """
    diz = diz["matrix"]
    sub_matrix = diz[n[:, None], m]
    sub_matrix = np.where(sub_matrix != np.inf, sub_matrix, 0)
    min_columns = np.amin(sub_matrix, axis=0)
    sum_columns = np.sum(min_columns)
    return sum_columns


def comparison_random_walk(network: nx.Graph, genesetA: list, genesetB: list, diz: dict = {}) -> float:
    """
    Evaluate the random walk on two genesets

    :param network: the graph representing the network
    :param genesetA: the first geneset list
    :param genesetB: the second geneset list
    :param diz: the dictionary containing the nodes name and index
    """
    try:
        diz["matrix"]
    except KeyError:
        print("The dictionary doesnt have a matrix key")
        raise

    genesetA_index = [diz["nodes"].index(i) for i in genesetA]
    genesetB_index = [diz["nodes"].index(i) for i in genesetB]

    if len(genesetA_index) == 0 or len(genesetB_index) == 0:
        sys.exit()

    probAB = [diz["matrix"][i, j] for i in genesetA_index for j in genesetB_index]
    probBA = [diz["matrix"][i, j] for i in genesetB_index for j in genesetA_index]

    prob = np.sum(probAB) + np.sum(probBA)
    return prob
