import abc

import pandas as pd
import tables
import networkx as nx
import numpy as np
import logging
import sys


class ReadingData(object):
    """Abstract class used to read different types of file. You can implement your own reading method, but remember
    that each subclass must implement the 'readfile' and 'get_data' methods
    """

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

    def __init__(self, filename: str, pd_table: bool = False, int_type: int = None):
        """
        :param filename: represents the path to the network file
        :param pd_table: if the results is going to be a pd.dataframe
        """
        super().__init__()
        self.filename = filename
        self.int_type = int_type
        self.pd_table = pd_table
        self.interactions = []

        if not self.pd_table:
            self.__readfile()
        self.graph = self._convert_to_graph()

    def __readfile(self) -> None:
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

    def _convert_to_graph(self) -> nx.Graph:
        """
        Converts the interactions into a graph object

        :return: graph from the interactions
        """
        graph = nx.Graph()
        graph.add_edges_from(self.interactions)
        graph.remove_edges_from(graph.selfloop_edges())
        return graph

    def get_data(self) -> pd.DataFrame or list:
        """
        Returns the data of the tsv file

        :return: list representing the genes read in the file
        """
        if self.pd_table:
            return pd.read_table(self.filename)
        else:
            return self.interactions

    def get_network(self) -> nx.Graph:
        """
        Returns the nx.graph object of the network

        :return: graph containing the network information
        """
        return self.graph


class ReadGmt(ReadingData):
    """
    A class used to read the .gmt gene file inside pygna
    """

    def __init__(self, filename: str, read_descriptor: bool = False):
        """
        :param filename: represents the path to the geneset file
        :param read_descriptor: if the descriptor is given. Default = False
        """
        super().__init__()
        self.filename = filename
        self.read_descriptor = read_descriptor

        self.gmt_data = self.__readfile()

    def __readfile(self) -> dict:
        """
        This method reads the geneset file into a variable

        :return: gene_list representing the list of genes
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

    def get_data(self) -> dict:
        """
        Returns the data of the gmt file

        :return: dict representing the genes list
        """
        return self.gmt_data

    def get_geneset(self, setname: str = None) -> dict:
        """
        Returns the geneset from the gmt file

        :param setname: str, the setname to extract
        :return: the geneset data
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

    def __init__(self, filename: str, sep: str = ",", use_cols: list = None, column_to_fill: str = None):
        """
        :param filename: represents the path to the data file
        :param sep: the separator to be used
        :param use_cols: columns used to be read and grouped
        :param column_to_fill: column to fill the NA values
        """
        super().__init__()
        self.filename = filename
        self.sep = sep
        self.use_cols = use_cols
        self.name_column = column_to_fill

        self.data = self.__readfile()
        if self.name_column is not None:
            self._fill_na_column()

    def __readfile(self) -> pd.DataFrame:
        """
        This method read the file and saves the data inside a class attribute

        :return: dataframe representing teh data read inside the .csv
        """
        with open(self.filename, "r") as f:
            table = pd.read_csv(f, sep=self.sep, usecols=self.use_cols)
            return table

    def get_data(self) -> pd.DataFrame:
        """
        Returns the data of the csv file

        :return: dataframe representing the data read inside the .csv
        """
        return self.data

    def _fill_na_column(self) -> None:
        """
        Fill the N/A values with a (str) 0
        """
        self.data[self.name_column].fillna(0, inplace=True)
        self.data[self.name_column] = self.data[self.name_column].astype(int)


class ReadTxt(ReadingData):
    """
    This class reads a txt file containing a single gene per line
    """

    def __init__(self, filename: str):
        super().__init__()
        self.filename = filename
        self.data = []

    def __readfile(self) -> None:
        """
        Read the file, line per line
        """
        with open(self.filename, "r") as f:
            gene_line = f.readline()
            while gene_line:
                self.data.append(gene_line)
                gene_line = f.readline()

    def get_data(self) -> pd.DataFrame:
        """
        Get the dataframe from the class

        :return: dataframe object from the file read
        """
        return pd.DataFrame(self.data)


class ReadDistanceMatrix(ReadingData):
    """
    This class read a distance matrix in the HDF5 format
    """

    def __init__(self, filename: str, in_memory: bool = False):
        """
        :param filename: the path of the file to be read
        :param in_memory: keep the matrix in memory or not
        """
        super().__init__()
        self.filename = filename
        self.nodes = None
        self.data = None
        self.memory = in_memory

        self.__readfile()
        if type(self.nodes[0]) == bytes:
            self._decode()

    def __readfile(self) -> None:
        """
        This method read and stores matrix information in memory or by reading it on the disk
        """
        if self.memory:
            hdf5_file = tables.open_file(self.filename, mode="r", driver="H5FD_CORE")
        else:
            hdf5_file = tables.open_file(self.filename, mode="r")
        self.nodes = list(hdf5_file.root.nodes[:])
        self.data = hdf5_file.root.matrix[:]

    def _decode(self) -> None:
        """
        Elaborate teh nodes from the graph
        """
        self.nodes = [i.decode() for i in self.nodes]

    def get_data(self) -> [list, np.matrix]:
        """Return the data of the HDF5 Matrix

        :return: table data, the data of the HDF5 Matrix and table nodes, the nodes of the HDF5 Matrix
        """
        return self.nodes, self.data
