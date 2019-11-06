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


def map_table2matrix(node_list, x):
    return node_list.index(x)


class DiffusionTest:
    def __init__(
        self,
        test_statistic,
        nodes,
        diffusion_matrix,
        weights_table,
        names_col="name",
        weights_col="stat",
        diz={},
    ):

        self.test_statistic = test_statistic
        self.nodes = nodes
        self.diffusion_matrix = diffusion_matrix

        self.diz = diz

        # print(type(self.__network))
        self.universe = set(self.nodes)

        self.table = weights_table[weights_table[names_col].isin(self.nodes)]
        logging.info("weights table: %d rows" % len(weights_table))
        logging.info("mapped table: %d mapped rows" % len(self.table))

        if len(self.table) < 10:
            logging.warning(
                "there are less than 10 elements in the table that are mapped to the universe"
            )

        self.table_index = [
            map_table2matrix(self.nodes, i)
            for i in self.table[names_col].values.tolist()
        ]

        self.weights = np.zeros((self.diffusion_matrix.shape[1], 1))

        print(self.weights.shape)
        for k in range(len(self.table_index)):
            self.weights[self.table_index[k], 0] = self.table[
                weights_col
            ].values.tolist()[k]

    def empirical_pvalue(self, geneset, alternative="less", max_iter=100, cores=1):
        # mapping geneset

        geneset = geneset
        mapped_geneset = list(self.universe.intersection(set(geneset)))

        if len(mapped_geneset) == 0:
            return 0, 0, np.array([0]), 0, 0
        else:
            geneset_index = [map_table2matrix(self.nodes, i) for i in mapped_geneset]
            logging.info(
                "Mapped %d genes out of %d." % (len(mapped_geneset), len(geneset))
            )

            observed = self.test_statistic(
                self.diffusion_matrix,
                self.weights,
                geneset_index,
                self.diz,
                observed_flag=True,
            )
            logging.info("Observed %f." % observed)

            # iterations
            null_distribution = DiffusionTest.get_null_distribution_mp(
                self, geneset_index, max_iter, n_proc=cores
            )
            # computing empirical pvalue
            pvalue = 1
            if alternative == "greater":
                pvalue = np.sum(null_distribution >= observed) / float(max_iter)
            else:
                pvalue = np.sum(null_distribution <= observed) / float(max_iter)

            return (
                observed,
                pvalue,
                null_distribution,
                len(mapped_geneset),
                len(geneset),
            )

    def get_null_distribution_mp(self, geneset_index, iter=100, n_proc=1):

        print("n_proc=" + str(n_proc))

        if n_proc == 1:
            null_distribution = DiffusionTest.get_null_distribution(
                self, geneset_index, iter
            )

        else:

            p = multiprocessing.Pool(n_proc)
            n_trial = int(iter / n_proc)
            print("n_trial=" + str(n_trial))
            results = [
                p.apply_async(
                    DiffusionTest.get_null_distribution,
                    args=(self, geneset_index, n_trial),
                )
                for w in list(range(1, n_proc + 1))
            ]
            null_distribution = np.array([])
            for r in results:
                null_distribution = np.hstack((null_distribution, np.array(r.get())))
            print(len(null_distribution))
            p.close()

        return np.asarray(null_distribution)

    def get_null_distribution(self, geneset_index, n_samples, randomize="index"):
        np.random.seed()
        random_dist = []
        for i in range(n_samples):
            random_weights = self.weights.copy()
            if randomize != "index":
                np.random.shuffle(random_weights)
            else:
                np.random.shuffle(geneset_index)
            random_dist.append(self.test_statistic(
                self.diffusion_matrix,
                random_weights,
                geneset_index,
                self.diz,
                observed_flag=True,
                )
            )
        return random_dist


###############################################################################
###  TEST STATISTICS FOR SINGLE GENESET  ######################################
###############################################################################


def weights_diffusion_statistic(
    matrix, weights, geneset_index, diz={}, observed_flag=False
):

    """

    """

    if matrix.shape[1] != weights.shape[0]:
        logging.warning("pass the right shape for the weights")
    try:
        product = np.matmul(matrix, weights)
    except:
        logging.warning("error in matrix multiplication")

    w = [product[i] for i in geneset_index]
    stat = np.sum(w) / len(w)
    return stat

def hotnet_diffusion_statistic(
    matrix, weights, geneset_index, diz={}, observed_flag=False
):

    """
    Applies the diagonal matrix of weights and gets all rows and
    columns according to the genelist
    """

    weights =np.diagflat(weights.T)
    if matrix.shape[1] != weights.shape[0]:
        logging.warning("pass the right shape for the weights")
    try:
        product = np.matmul(matrix, weights)
    except:
        logging.warning("error in matrix multiplication")

    prob = [product[i,j] for i in geneset_index for j in geneset_index if i != j]
    stat = np.sum(prob)  
    return stat