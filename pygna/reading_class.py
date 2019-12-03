import abc

import pandas as pd
import tables
import networkx as nx
import logging
import sys


class ReadingData(object):
    """Abstract class used to read different types of file. Each subclass must implement the 'readfile'
    and get_data method"""
    def __init__(self):
        super(ReadingData, self).__init__()

    @abc.abstractmethod
    def __readfile(self):
        raise NotImplementedError

    @abc.abstractmethod
    def get_data(self):
        raise NotImplementedError


class ReadTsv(ReadingData):
    """
    A class used to read the .tsv network file inside pygna
    """

    def __init__(self, filename, int_type=None):
        # TODO Fix documentation about int_type
        """
        :param filename: str, represents the path to the network file
        :param int_type: Unknown
        """
        super().__init__()
        self.filename = filename
        self.int_type = int_type
        self.interactions = []

        self.__readfile()
        self.graph = self._convert_to_graph()

    def __readfile(self):
        """
        This method read the file and saves the data inside a class attribute
        """
        with open(self.filename, "r") as f:
            for record in f:
                if record.startswith("#"):
                    continue

                fields = record.strip().split("\t")
                if self.int_type:
                    types = fields[3].split(";")
                    if self.int_type in types:
                        self.interactions.append((fields[0], fields[1]))
                    else:
                        continue
                else:
                    self.interactions.append((fields[0], fields[1]))

    def _convert_to_graph(self):
        """
        Converts the interactions into a graph object
        :return: nx.graph, from the interactions
        """
        graph = nx.Graph()
        graph.add_edges_from(self.interactions)
        graph.remove_edges_from(graph.selfloop_edges())
        return graph

    def get_data(self):
        """
        Returns the data of the tsv file
        :return: list, represents the genes read in the file
        """
        return self.interactions

    def get_pd_data(self):
        """
        Returns the data into a pandas object
        :return: pd.dataframe, represents the data into a pandas object
        """
        return pd.read_table(self.filename)

    def get_network(self):
        """
        Returns the nx.graph object of the network
        :return: nx.graph, containing the network information
        """
        return self.graph


class ReadGmt(ReadingData):
    """
    A class used to read the .gmt gene file inside pygna
    """

    def __init__(self, filename, read_descriptor=False):
        """
        :param filename: str, represents the path to the geneset file
        :param read_descriptor: bool, if the descriptor is given. Default = False
        """
        super().__init__()
        self.filename = filename
        self.read_descriptor = read_descriptor

        self.gmt_data = self.__readfile()

    def __readfile(self):
        """
        This method reads the geneset file into a variable
        :returns gene_list: dict, represents the genes list
        """
        gene_lists = dict()
        with open(self.filename, "r") as f:
            for record in f:
                fields = record.strip().split("\t")
                if self.read_descriptor:
                    gene_lists[fields[0]] = {}
                    gene_lists[fields[0]]["genes"] = fields[2:]
                    gene_lists[fields[0]]["descriptor"] = fields[1]
                else:
                    gene_lists[fields[0]] = fields[2:]
            return gene_lists

    def get_data(self):
        """
        Returns the data of the gmt file
        :return: dict, represents the genes list
        """
        return self.gmt_data

    def get_geneset(self, setname):
        """
        Returns the geneset from the gmt file
        :param setname: str, the setname to extract
        :return: pd.dataframe, the geneset data
        """
        if setname is not None:
            if setname in self.gmt_data:
                temp = self.gmt_data[setname]
                self.gmt_data.clear()
                self.gmt_data[setname] = temp
            else:
                logging.error("Cannot find geneset: %s" % setname)
                sys.exit(-1)
        return self.gmt_data


class ReadCsv(ReadingData):
    """
    A class used to read the .csv data file
    """

    def __init__(self, filename, sep=",", use_cols=None):
        """
        :param filename: str, represents the path to the data file
        :param sep: str, the separator to be used
        :param use_cols: list, columns used to be read and grouped
        """
        super().__init__()
        self.filename = filename
        self.sep = sep
        self.use_cols = use_cols

        self.data = self.__readfile()

    def __readfile(self):
        """
        This method read the file and saves the data inside a class attribute
        :return: pd.dataframe, represents teh data read inside the .csv
        """
        with open(self.filename, "r") as f:
            table = pd.read_csv(f, sep=self.sep, usecols=self.use_cols)
            return table

    def get_data(self):
        """
        Returns the data of the csv file
        :return: pd.dataframe, represents teh data read inside the .csv
        """
        return self.data

    def fill_na_column(self, name_column="gene_name"):
        """
        Fill the N/A values with a (str) 0
        :param name_column: name of the column to filter
        :return null
        """
        self.data[name_column] = self.data[name_column].fillna(0).apply(int).apply(str)


class ReadDistanceMatrix(ReadingData):
    """
    This class read a distance matrix in the HDF5 format
    """
    def __init__(self, filename, in_memory=False):
        """
        :param filename: str, the path of the file to be read
        """
        super().__init__()
        self.filename = filename
        self.nodes = None
        self.data = None
        self.memory = in_memory

        self.__readfile()
        if type(self.nodes[0]) == bytes:
            self._decode()

    def __readfile(self):
        """
        This method read and stores matrix information in memory or by reading it on the disk
        """
        if self.memory:
            hdf5_file = tables.open_file(self.filename, mode="r", driver="H5FD_CORE")
        else:
            hdf5_file = tables.open_file(self.filename, mode="r")
        self.nodes = list(hdf5_file.root.nodes[:])
        self.data = hdf5_file.root.matrix[:]

    def _decode(self):
        """
        Elaborate teh nodes from the graph
        :return: null
        """
        self.nodes = [i.decode() for i in self.nodes]

    def get_data(self):
        """
        Return the data of the HDF5 Matrix
        :return: table data, the data of the HDF5 Matrix
        :return: table nodes, the nodes of the HDF5 Matrix
        """
        return self.nodes, self.data
