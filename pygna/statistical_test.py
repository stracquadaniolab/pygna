import logging
import networkx as nx
import numpy as np
import sys
import multiprocessing


class StatisticalTest:
    """
    This class implements the statistical analysis performed by Pygna.
    It performs the statistical tests on the given network, elaborates the number of observed genes, the pvalue etc.
    Please refer to the single method documentation for the returning values
    """
    def __init__(self, test_statistic, network: nx.Graph, diz: dict = {}):
        """
        :param test_statistic: the statistical function to be used for the calculation of the empirical p-value and the
        null distribution
        :param network: the network to be used for the analysis
        :param diz: the dictionary containing the genes
        """
        self.__test_statistic = test_statistic
        self.__network = network
        self.__diz = diz

        if (type(self.__network) is nx.Graph) or (type(self.__network) is nx.DiGraph):
            self.__universe = set(self.__network.nodes())
        elif type(self.__network) is dict:
            self.__universe = set(self.__network.keys())
        else:
            logging.error("Unknown network type: %s" % type(self.__network))
            sys.exit(-1)

    def empirical_pvalue(self, geneset: set, alternative: str = "less", max_iter: int = 100, cores: int = 1) -> \
        [int, float, float, int, int]:
        """
        Calculate the empirical pvalue on the genes list

        :param geneset: the geneset to elaborate
        :param alternative: the pvalue selection of the observed genes
        :param max_iter: the number of iterations to be performed
        :param cores: the number of cores to be used
        :return observed, pvalue, null_distribution, len(mapped_genesetA), len(mapped_genesetB): the list with the data calculated
        """
        # mapping geneset
        mapped_geneset = sorted(list(set(geneset).intersection(self.__universe)))
        if len(mapped_geneset) == 0:
            return 0, 0, np.array([0]), 0, 0
        else:
            logging.info("Mapped %d genes out of %d." % (len(mapped_geneset), len(geneset)))
            observed = self.__test_statistic(self.__network, mapped_geneset, self.__diz, observed_flag=True)
            # iterations
            null_distribution = StatisticalTest.get_null_distribution_mp(self, mapped_geneset, max_iter, n_proc=cores)
            # computing empirical pvalue
            pvalue = 1
            if alternative == "greater":
                pvalue = (np.sum(null_distribution >= observed) + 1) / (float(max_iter) + 1)
            else:
                pvalue = (np.sum(null_distribution <= observed) + 1) / (float(max_iter) + 1)

            return (
                observed,
                pvalue,
                null_distribution,
                len(mapped_geneset),
                len(geneset),
            )

    def get_null_distribution_mp(self, geneset: list, iter: int = 100, n_proc: int = 1):
        """
        Calculate the null distribution with multiple cores on the geneset

        :param geneset: the geneset to be used
        :param iter: the number of iterations to perform
        :param n_proc: the number of cpu to use for the elaboration
        :return: the array with null distribution
        """
        if n_proc == 1:
            null_distribution = StatisticalTest.get_null_distribution(self, geneset, iter)

        else:

            p = multiprocessing.Pool(n_proc)
            n_trial = int(iter / n_proc)
            results = [p.apply_async(StatisticalTest.get_null_distribution,
                                     args=(self, geneset, n_trial)) for w in list(range(1, n_proc + 1))]
            null_distribution = np.array([])
            for r in results:
                null_distribution = np.hstack((null_distribution, np.array(r.get())))
            p.close()

        return np.asarray(null_distribution)

    def get_null_distribution(self, geneset: list, n_samples: int):
        """
        Calculate the null distribution over the geneset

        :param geneset: the geneset to be used
        :param n_samples: the number of samples to be taken
        :return: the random distribution calculated
        """
        np.random.seed()
        random_dist = []
        for i in range(n_samples):
            random_sample = np.random.choice(list(self.__universe), len(geneset), replace=False)
            random_dist.append(self.__test_statistic(self.__network, set(random_sample), self.__diz))

        return random_dist


###############################################################################
###  TEST STATISTICS FOR SINGLE GENESET  ######################################
###############################################################################


def geneset_localisation_statistic_median(network: nx.Graph, geneset: set, diz: dict = {}, observed_flag: bool = False) \
    -> float:
    """
    Calculate the median shortest path for each node

    :param network: the network used in the analysis
    :param geneset: the geneset to analyse
    :param diz: the dictionary containing the genes name
    :param observed_flag: XXXX
    """
    cum_sum = 0.0
    geneset_index = [diz["nodes"].index(i) for i in geneset]

    for u in geneset_index:
        d_uv = []
        for v in geneset_index:
            d_uv.append(diz["matrix"][v][u] + diz["matrix"][u][v])

        cum_sum += np.median(d_uv)
    return cum_sum / float(len(geneset))


def geneset_localisation_statistic(network: nx.Graph, geneset: set, diz: dict = {}, observed_flag: bool = False) \
    -> float:
    """
    Identify the genes in a geneset

    :param network: the network used in the analysis
    :param geneset: the geneset to analyse
    :param diz: the dictionary containing the genes name
    :param observed_flag: XXXX
    """
    n = np.array([diz["nodes"].index(i) for i in geneset])
    diz = diz["matrix"]

    sub_matrix = diz[n[:, None], n]
    min_columns = np.amin(sub_matrix, axis=0)
    sum_columns = np.sum(min_columns)
    return sum_columns / len(n)


def geneset_module_statistic(network: nx.Graph, geneset: set, diz: dict = {}, observed_flag: bool = False) -> float:
    """
    Evaluate the length of a observed network

    :param network: the network used in the analysis
    :param geneset: the geneset to analyse
    :param diz: the dictionary containing the genes name
    :param observed_flag: XXXX
    """
    # Largest Connected Component for the subgraph induced by the geneset
    module = nx.subgraph(network, geneset)

    if observed_flag == True:
        pass
    cc = sorted(list(nx.connected_components(module)), key=len, reverse=True)

    if len(cc) > 0:
        return len(cc[0])
    else:
        return 0


def geneset_total_degree_statistic(network: nx.Graph, geneset: set, diz: dict = {}, observed_flag: bool = False) -> \
    float:
    """
    Total degree of the geneset: average total_degree

    :param network: the network used in the analysis
    :param geneset: the geneset to analyse
    :param diz: the dictionary containing the genes name
    :param observed_flag: XXXX
    """
    degree = nx.degree(network)
    geneset = list(geneset)
    total = np.array([degree[g] for g in geneset])
    return np.average(total)


def geneset_internal_degree_statistic(network: nx.Graph, geneset: set, diz: dict = {}, observed_flag: bool = False) \
    -> float:
    """
    Internal degree ratio: average of the ratio internal_degree/total_degree

    :param network: the network used in the analysis
    :param geneset: the geneset to analyse
    :param diz: the dictionary containing the genes name
    :param observed_flag: XXXX
    """
    degree = nx.degree(network)
    total = np.array([degree[g] for g in geneset])

    subgraph = network.subgraph(geneset)
    degree_internal = nx.degree(subgraph)
    internal = np.array([degree_internal[g] for g in geneset])

    ratio = internal / total
    ratio[total == 0] = 0.5

    return np.average(ratio)


def geneset_RW_statistic(network: nx.Graph, geneset: set, diz: dict = {}, observed_flag: bool = False) -> np.ndarray:
    """
    Poisson binomial probability, sum of interaction probabilities for the genes in the geneset

    :param network: the network used in the analysis
    :param geneset: the geneset to analyse
    :param diz: the dictionary containing the genes name
    :param observed_flag: XXXX
    """
    try:
        diz["matrix"]
    except KeyError:
        print("The dictionary doesnt have a matrix key")
        raise

    geneset_index = [diz["nodes"].index(i) for i in geneset]
    prob = [diz["matrix"][i, j] for i in geneset_index for j in geneset_index if i != j]
    prob = np.sum(prob)
    return prob
