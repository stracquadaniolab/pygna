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

    def __init__(self, filename, pd_table=False, int_type=None):
        # TODO Fix documentation about int_type
        """
        :param filename: str, represents the path to the network file
        :param pd_table: bool, if the results is going to be a pd.dataframe
        :param int_type: Unknown
        """
        super().__init__()
        self.filename = filename
        self.int_type = int_type
        self.pd_table = pd_table
        self.interactions = []

        if not self.pd_table:
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
        if self.pd_table:
            return pd.read_table(self.filename)
        else:
            return self.interactions

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

    def get_geneset(self, setname=None):
        """
        Returns the geneset from the gmt file
        :param setname: str, the setname to extract
        :return: pd.dataframe, the geneset data
        """
        # TODO The following "if" is because the parameter "setname" is defined in the calling functions above this one.
        #  This should not be allowed to keep a good encapsulation of this method. After the refactor, this should be
        #  check
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

    def __init__(self, filename, sep=",", use_cols=None, column_to_fill=None):
        """
        :param filename: str, represents the path to the data file
        :param sep: str, the separator to be used
        :param use_cols: list, columns used to be read and grouped
        :param column_to_fill: str, column to fill the NA values
        """
        super().__init__()
        self.filename = filename
        self.sep = sep
        self.use_cols = use_cols
        self.name_column = column_to_fill

        self.data = self.__readfile()
        if self.name_column is not None:
            self._fill_na_column()

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

    def _fill_na_column(self):
        """
        Fill the N/A values with a (str) 0
        :return null
        """
        self.data[self.name_column].fillna(0, inplace=True)
        self.data[self.name_column] = self.data[self.name_column].astype(int)


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


# class ReadNnDataset(ReadingData):
#     def __init__(self, filename, dataset_name):
#         super(ReadNnDataset, self).__init__()
#         self.filename = filename
#         self.dataset_name = dataset_name
#
#         self.adj, self.features, self.labels, self.idx_train, self.idx_val, self.idx_test = self.__readfile()
#
#     def __readfile(self):
#         """Load citation network dataset (cora only for now)"""
#         logging.info('Loading {} dataset...'.format(self.dataset_name))
#
#         idx_features_labels = np.genfromtxt("{}{}.content".format(self.filename, self.dataset_name),
#                                             dtype=np.dtype(str))
#         features = sp.csr_matrix(idx_features_labels[:, 1:-1], dtype=np.float32)
#         labels = self._encode_onehot(idx_features_labels[:, -1])
#
#         # build graph
#         idx = np.array(idx_features_labels[:, 0], dtype=np.int32)
#         idx_map = {j: i for i, j in enumerate(idx)}
#         edges_unordered = np.genfromtxt("{}{}.cites".format(self.filename, self.dataset_name), dtype=np.int32)
#         edges = np.array(list(map(idx_map.get, edges_unordered.flatten())), dtype=np.int32).reshape(
#             edges_unordered.shape)
#         adj = sp.coo_matrix((np.ones(edges.shape[0]), (edges[:, 0], edges[:, 1])),
#                             shape=(labels.shape[0], labels.shape[0]), dtype=np.float32)
#
#         # build symmetric adjacency matrix
#         adj = adj + adj.T.multiply(adj.T > adj) - adj.multiply(adj.T > adj)
#
#         features = self._normalize_features(features)
#         adj = self._normalize_adj(adj + sp.eye(adj.shape[0]))
#
#         idx_train = range(140)
#         idx_val = range(200, 500)
#         idx_test = range(500, 1500)
#
#         adj = torch.FloatTensor(np.array(adj.todense()))
#         features = torch.FloatTensor(np.array(features.todense()))
#         labels = torch.LongTensor(np.where(labels)[1])
#
#         idx_train = torch.LongTensor(idx_train)
#         idx_val = torch.LongTensor(idx_val)
#         idx_test = torch.LongTensor(idx_test)
#
#         return adj, features, labels, idx_train, idx_val, idx_test
#
#     def _normalize_adj(self, mx):
#         """Row-normalize sparse matrix"""
#         rowsum = np.array(mx.sum(1))
#         r_inv_sqrt = np.power(rowsum, -0.5).flatten()
#         r_inv_sqrt[np.isinf(r_inv_sqrt)] = 0.
#         r_mat_inv_sqrt = sp.diags(r_inv_sqrt)
#         return mx.dot(r_mat_inv_sqrt).transpose().dot(r_mat_inv_sqrt)
#
#     def _normalize_features(self, mx):
#         """Row-normalize sparse matrix"""
#         rowsum = np.array(mx.sum(1))
#         r_inv = np.power(rowsum, -1).flatten()
#         r_inv[np.isinf(r_inv)] = 0.
#         r_mat_inv = sp.diags(r_inv)
#         mx = r_mat_inv.dot(mx)
#         return mx
#
#     def _encode_onehot(self, labels):
#         classes = set(labels)
#         classes_dict = {c: np.identity(len(classes))[i, :] for i, c in enumerate(classes)}
#         labels_one_hot = np.array(list(map(classes_dict.get, labels)), dtype=np.int32)
#         return labels_one_hot
#
#     def get_data(self):
#         return self.adj, self.features, self.labels, self.idx_train, self.idx_val, self.idx_test
