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


class StatisticalTest:
    def __init__(self, test_statistic, network, diz={}):
        self.__test_statistic = test_statistic

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

    def empirical_pvalue(self, geneset, alternative="less", max_iter=100, cores=1):
        # mapping geneset
        mapped_geneset = sorted(list(set(geneset).intersection(self.__universe)))
        if len(mapped_geneset) == 0:
            return 0, 0, np.array([0]), 0, 0
        else:
            logging.info(
                "Mapped %d genes out of %d." % (len(mapped_geneset), len(geneset))
            )
            observed = self.__test_statistic(
                self.__network, mapped_geneset, self.__diz, observed_flag=True
            )
            # TODO: remove begin
            # verbose_localisation_statistic(self.__network, mapped_geneset, self.__diz, observed_flag=True)
            # TODO: remove end
            logging.info("Observed %f." % observed)
            # iterations
            null_distribution = StatisticalTest.get_null_distribution_mp(
                self, mapped_geneset, max_iter, n_proc=cores
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

    def get_null_distribution_mp(self, geneset, iter=100, n_proc=1):

        print("n_proc=" + str(n_proc))

        if n_proc == 1:
            null_distribution = StatisticalTest.get_null_distribution(
                self, geneset, iter
            )

        else:

            p = multiprocessing.Pool(n_proc)
            n_trial = int(iter / n_proc)
            print("n_trial=" + str(n_trial))
            results = [
                p.apply_async(
                    StatisticalTest.get_null_distribution, args=(self, geneset, n_trial)
                )
                for w in list(range(1, n_proc + 1))
            ]
            null_distribution = np.array([])
            for r in results:
                null_distribution = np.hstack((null_distribution, np.array(r.get())))
            print(len(null_distribution))
            p.close()

        return np.asarray(null_distribution)

    def get_null_distribution(self, geneset, n_samples):

        np.random.seed()
        random_dist = []
        for i in range(n_samples):
            random_sample = np.random.choice(
                list(self.__universe), len(geneset), replace=False
            )
            random_dist.append(
                self.__test_statistic(self.__network, set(random_sample), self.__diz)
            )

        return random_dist


###############################################################################
###  TEST STATISTICS FOR SINGLE GENESET  ######################################
###############################################################################


def geneset_localisation_statistic_median(
    network, geneset, diz={}, observed_flag=False
):
    """ median shortest path for each node """
    cum_sum = 0.0
    geneset_index = [diz["nodes"].index(i) for i in geneset]

    for u in geneset_index:
        d_uv = []
        for v in geneset_index:
            # print(v)
            d_uv.append(diz["matrix"][v][u] + diz["matrix"][u][v])

        cum_sum += np.median(d_uv)

    # logging.info("localisation statistic= %f", (cum_sum/float(len(geneset))))
    return cum_sum / float(len(geneset))


def geneset_localisation_statistic(network, geneset, diz={}, observed_flag=False):
    # minumum shortest path
    cum_sum = 0.0
    geneset_index = [diz["nodes"].index(i) for i in geneset]

    for u in geneset_index:
        min_du = float("inf")
        for v in geneset_index:
            # print(v)
            d_uv = diz["matrix"][v][u] + diz["matrix"][u][v]
            if u != v and d_uv < min_du:
                min_du = d_uv
        cum_sum += min_du

    # logging.info("localisation statistic= %f", (cum_sum/float(len(geneset))))
    return cum_sum / float(len(geneset))


def geneset_module_statistic(network, geneset, diz={}, observed_flag=False):

    # Largest Connected Component for the subgraph induced by the geneset
    module = nx.subgraph(network, geneset)

    if observed_flag == True:
        pass
        # nx.write_graphml(module, "../results/graph.graphml")
    cc = sorted(list(nx.connected_components(module)), key=len, reverse=True)
    # if observed_flag==True:
    #    diag.connected_components_diagnostic(cc)

    if len(cc) > 0:
        return len(cc[0])
    else:
        return 0


def geneset_total_degree_statistic(network, geneset, diz={}, observed_flag=False):
    """
    Total degree of the geneset: average total_degree
    """
    degree = nx.degree(network)
    geneset = list(geneset)
    total = np.array([degree[g] for g in geneset])
    # print(total)
    return np.average(total)


def geneset_internal_degree_statistic(network, geneset, diz={}, observed_flag=False):
    """
    Internal degree ratio: average of the ratio internal_degree/total_degree
    """
    degree = nx.degree(network)
    total = np.array([degree[g] for g in geneset])
    # print(total)

    subgraph = network.subgraph(geneset)
    degree_internal = nx.degree(subgraph)
    internal = np.array([degree_internal[g] for g in geneset])

    ##print(degree_internal)

    ratio = internal / total
    ratio[total == 0] = 0.5

    return np.average(ratio)


def geneset_RW_statistic(network, geneset, diz={}, observed_flag=False):

    """ Poisson binomial probability, sum of interaction probabilities for the genes in the geneset
    """

    try:
        diz["matrix"]
    except KeyError:
        print("The dictionary doesnt have a matrix key")
        raise

    geneset_index = [diz["nodes"].index(i) for i in geneset]
    prob = [diz["matrix"][i, j] for i in geneset_index for j in geneset_index if i != j]
    prob = np.sum(prob)  # /len(prob)
    return prob


##############################################################
################# UNUSED #####################################
##############################################################


def geneset_hitting_time_statistic(network, geneset, diz={}, observed_flag=False):

    try:
        diz["hitting_time"]
    except KeyError:
        print(
            "The hitting time dictionary must be passed as the value of a dictionary with key hitting_time"
        )
        raise

    nodes = diz["hitting_time"]["nodes"]

    hitting_time = [
        diz["hitting_time"]["tau"][nodes.index(i), nodes.index(j)]
        for i in geneset
        for j in geneset
        if i != j
    ]
    average_ht = np.sum(hitting_time) / len(hitting_time)
    # print("average")
    # print("Hitting time calculated")
    # print(average_ht)
    return average_ht


def geneset_neighbor_statistic(network, geneset, diz={}, observed_flag=False):
    # number of neighbors
    sg = set(geneset)
    num_edges = 0
    for g in geneset:
        num_edges += len(
            set(network.neighbors(g))
        )  # len(set(network.neighbors(g)) & sg)
    return num_edges


def verbose_localisation_statistic(network, geneset, diz={}, observed_flag=False):
    # minumum shortest path
    cum_sum = 0.0
    geneset_index = [diz["nodes"].index(i) for i in geneset]
    print(geneset_index)

    gu = -1
    for u in geneset_index:
        gu = gu + 1
        min_du = float("inf")
        closest = ""
        gv = -1
        for v in geneset_index:
            gv = gv + 1

            kkk = nx.shortest_path_length(
                network, source=geneset[gu], target=geneset[gv]
            )
            d_uv = diz["matrix"][v][u] + diz["matrix"][u][v]
            if kkk != d_uv:
                print("WARNING %d not equal %d" % (kkk, d_uv))

            if u != v and d_uv < min_du:
                min_du = d_uv
                closest = v

        print("For gene %s min_du is: %d , gene %s" % (str(u), min_du, str(closest)))
        cum_sum += min_du

    # logging.info("localisation statistic= %f", (cum_sum/float(len(geneset))))
    return cum_sum / float(len(geneset))

