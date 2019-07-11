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
