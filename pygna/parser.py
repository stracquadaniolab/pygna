import logging
import pickle
import networkx as nx
import tables
import pandas as pd

class Parser:

    def __init__(self):
        pass

    def read(self, filename, organism='9606', int_type = None):
        pass


class Tab2Parser(Parser):

    __FIELD_ENTREZ_A__ = 7
    __FIELD_ENTREZ_B__ = 8
    __FIELD_ORGANISM_A__ = 15
    __FIELD_ORGANISM_B__ = 16

    def __init__(self):
        pass

    def read(self, filename, organism='9606', int_type = None):
        interactions = []
        for record in open(filename):
            fields = record.strip().split("\t")
            if fields[Tab2Parser.__FIELD_ORGANISM_A__] == organism and fields[Tab2Parser.__FIELD_ORGANISM_B__] == organism:
                gene_a = fields[Tab2Parser.__FIELD_ENTREZ_A__]
                gene_b = fields[Tab2Parser.__FIELD_ENTREZ_B__]
                interactions.append((gene_a, gene_b))

        graph = nx.Graph()
        graph.add_edges_from(interactions)
        graph.remove_edges_from(graph.selfloop_edges())
        return graph

class TSVParser(Parser):
    def read(self, filename, organism='9606', int_type = None):
        interactions = []
        for record in open(filename):
            if record.startswith('#'):
                continue
            
            fields = record.strip().split("\t")
            if int_type:
                types=fields[3].split(";")
                if int_type in types:
                    interactions.append((fields[0], fields[1]))                
                else:
                    continue
            else:
                interactions.append((fields[0], fields[1]))

        graph = nx.Graph()
        graph.add_edges_from(interactions)
        graph.remove_edges_from(graph.selfloop_edges())
        return graph



class GMTParser(Parser):

    def read(self, filename, read_descriptor=False, organism='9606', int_type = None):
        genelists = dict()
        if read_descriptor:
            for record in open(filename):
                fields = record.strip().split('\t')
                genelists[fields[0]]={}
                genelists[fields[0]]["genes"]= fields[2:]
                genelists[fields[0]]["descriptor"]= fields[1]

        else:
            for record in open(filename):
                fields = record.strip().split('\t')
                genelists[fields[0]] = fields[2:]

        return genelists


class DistanceMatrixParser(Parser):
    def read(self, filename, organism='9606', int_type = None):
        return pickle.load(open(filename, 'rb'))

class DiffusionDictParser(Parser):
    def read(self, filename):
        return pickle.load(open(filename, 'rb'))

class Hdf5MatrixParser(Parser):

    def read(self, filename):

        with tables.open_file(filename, mode="r") as hdf5_file:
            hdf5_nodes=list(hdf5_file.root.nodes[:])
            hdf5_data=hdf5_file.root.matrix[:]

        return (hdf5_nodes, hdf5_data)

