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

import pygna.parser as ps
import pygna.reading_class as rc
'''
class Data:
    def __init__(self, geneset_file, setname = None, geneset_b_file=None, setname_b=None):

        self.__geneset_file = geneset_file
        self.__setname = setname
        self.__geneset = ps.__load_geneset(geneset_file, setname)
        self.__geneset_b = None
        self.to_universe = False

        if geneset_b_file:
            self.__geneset_b = ps.__load_geneset(geneset_b_file, setname_b)

        if ((geneset_b_file == None) & (type(setname_b)==string)):
            self.__geneset_b = ps.__load_geneset(geneset_file, setname_b)


    def prepare_genesets(self, ids2idx_map, universe):
        universe = set(universe)
        for setname in self.__geneset.keys():
            genes = list(universe.intersection(set(self.__geneset[setname]['genes'])))
            self.__geneset[setname]['mapped'] = map(lambda x: (self.ids2int_map[x]), genes)

        if self.geneset_b:
            for setname in self.__geneset_b.keys():
                genes = list(universe.intersection(set(self.__geneset_b[setname]['genes'])))
                self.__geneset_b[setname]['mapped'] = map(lambda x: (self.ids2int_map[x]), genes)

        self.to_universe = True


class Network:
    def __init__(self, network_file, lcc = True, use_mapping = True):
        self.network
    network, ids2int_dict = rc.ReadTsv(network_file, use_mapping=True).get_network()
    network = nx.Graph(network.subgraph(max(nx.connected_components(network), key=len)))

'''
