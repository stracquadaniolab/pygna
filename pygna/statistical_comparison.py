import logging
import random
import networkx as nx
import numpy as np
import scipy
import sys
import matplotlib.pyplot as plt
import pygna.diagnostic as diag
import multiprocessing
import time
import pygna.statistical_test as st


class StatisticalComparison:
    def __init__(self, comparison_statistic, network, n_proc=1, diz={}):
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

    def comparison_empirical_pvalue(
        self, genesetA, genesetB, alternative="less", max_iter=100, keep=False
    ):

        """
        This method applies the __comparison_statistic to two genesets and
        returns a p-value
        """

        # mapping genesets
        mapped_genesetA = sorted(list(set(genesetA) & self.__universe))
        mapped_genesetB = sorted(list(set(genesetB) & self.__universe))
        logging.info(
            "Mapped %d genes out of %d from genesetA"
            % (len(mapped_genesetA), len(genesetA))
        )
        logging.info(
            "Mapped %d genes out of %d from genesetB"
            % (len(mapped_genesetB), len(genesetB))
        )

        # Observed statistical values
        observed = self.__comparison_statistic(
            self.__network, mapped_genesetA, mapped_genesetB, self.__diz
        )

        # iterations
        null_distribution = StatisticalComparison.get_comparison_null_distribution_mp(
            self, mapped_genesetA, mapped_genesetB, max_iter, keep
        )
        pvalue = 1
        if alternative == "greater":
            pvalue = np.sum(null_distribution >= observed) / float(
                len(null_distribution)
            )
        else:
            pvalue = np.sum(null_distribution <= observed) / float(
                len(null_distribution)
            )

        return (
            observed,
            pvalue,
            null_distribution,
            len(mapped_genesetA),
            len(mapped_genesetB),
        )

    def get_comparison_null_distribution_mp(
        self, genesetA, genesetB, max_iter=100, keep=False
    ):

        n_trial = int(max_iter / self.__n_proc)
        logging.info(
            "n_proc = %d, each computing %d permutations "
            % (int(self.__n_proc), n_trial)
        )
        
        # Don't set up the multicore architecture if only one core is given
        if self.__n_proc == 1:
            
            null_distribution = StatisticalComparison.get_comparison_null_distribution(
                self, genesetA, genesetB, max_iter, keep
            )

        else:
            
            p = multiprocessing.Pool(self.__n_proc)

            results = [
                p.apply_async(
                    StatisticalComparison.get_comparison_null_distribution,
                    args=(self, genesetA, genesetB, n_trial, keep),
                )
                for w in list(range(1, self.__n_proc + 1))
            ]

            null_distribution = np.array([])
            for r in results:
                null_distribution = np.hstack((null_distribution, np.array(r.get())))
            p.close()

        return np.asarray(null_distribution)

    def get_comparison_null_distribution(self, genesetA, genesetB, n_samples, keep):

        np.random.seed()
        random_dist = []

        # association statistic, B kept, A bootstrapped
        if keep:
            for i in range(n_samples):
                random_sample_A = np.random.choice(
                    list(self.__universe), len(genesetA), replace=False
                )
                random_dist.append(
                    self.__comparison_statistic(
                        self.__network, set(random_sample_A), set(genesetB), self.__diz
                    )
                )
        # association statistic, B kept, A bootstrapped
        else:
            for i in range(n_samples):
                random_sample_A = np.random.choice(
                    list(self.__universe), len(genesetA), replace=False
                )
                random_sample_B = np.random.choice(
                    list(self.__universe), len(genesetB), replace=False
                )
                random_dist.append(
                    self.__comparison_statistic(
                        self.__network,
                        set(random_sample_A),
                        set(random_sample_B),
                        self.__diz,
                    )
                )
        return random_dist


###############################################################################
###  TEST STATISTICS FOR COMPARISONS  #########################################
###############################################################################


def comparison_shortest_path(network, genesetA, genesetB, diz={}):

    cum_sum = 0.0

    genesetA_index = [diz["nodes"].index(i) for i in genesetA]
    genesetB_index = [diz["nodes"].index(i) for i in genesetB]

    if len(genesetA_index) == 0 or len(genesetB_index) == 0:
        sys.exit()

    for u in genesetA_index:
        min_du = float("inf")
        for v in genesetB_index:
            d_uv = diz["matrix"][v][u] + diz["matrix"][u][v]
            if u != v and d_uv < min_du:
                min_du = d_uv
            elif u == v:
                min_du = 0

        cum_sum += min_du

    cum_sum_v = 0
    for v in genesetB_index:
        min_dv = float("inf")
        for u in genesetA_index:
            d_uv = diz["matrix"][v][u] + diz["matrix"][u][v]
            if u != v and d_uv < min_dv:
                min_dv = d_uv
            elif u == v:
                min_dv = 0
        cum_sum_v += min_dv

        cum_sum += min_dv

    d_AB = cum_sum / float(len(genesetA) + len(genesetB))

    d_A = st.geneset_localisation_statistic(network, genesetA, diz)
    d_B = st.geneset_localisation_statistic(network, genesetB, diz)
    return d_AB - (d_A + d_B) / 2


def comparison_random_walk(network, genesetA, genesetB, diz={}):

    try:
        diz["matrix"]
    except KeyError:
        print("The dictionary doesnt have a matrix key")
        raise

    genesetA_index = [diz["nodes"].index(i) for i in genesetA]
    genesetB_index = [diz["nodes"].index(i) for i in genesetB]

    if len(genesetA_index) == 0 or len(genesetB_index) == 0:
        sys.exit()

    probAB = [
        diz["matrix"][i, j] for i in genesetA_index for j in genesetB_index
    ]  
    probBA = [
        diz["matrix"][i, j] for i in genesetB_index for j in genesetA_index
    ]  

    prob = np.sum(probAB) + np.sum(probBA)  
    return prob
