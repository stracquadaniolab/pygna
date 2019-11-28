from abc import ABC, abstractmethod

import pandas as pd


class ReadingData(ABC):
    """Abstract class used to read different types of file. Each subclass must implement the 'readfile'
    and get_data method"""
    def __init__(self):
        super(ReadingData, self).__init__()

    @abstractmethod
    def __readfile(self, filename):
        raise NotImplementedError

    @abstractmethod
    def get_data(self):
        raise NotImplementedError


class ReadTsv(ReadingData, ABC):
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
        super().__init__(filename)
        self.filename = filename
        self.int_type = int_type
        self.pd_table = pd_table

        if not self.pd_table:
            self.interactions = self._ReadingData__readfile(self.filename)

    def _ReadingData__readfile(self, filename):
        """
        This method read the file and saves the data inside a class attribute
        :param filename: str, represents the path to the file
        :return: interactions: list, represents the genes read in the file
        """
        interactions = []
        with open(filename, "r") as f:
            for record in f:
                if record.startswith("#"):
                    continue

                fields = record.strip().split("\t")
                if self.int_type:
                    types = fields[3].split(";")
                    if self.int_type in types:
                        interactions.append((fields[0], fields[1]))
                    else:
                        continue
                else:
                    interactions.append((fields[0], fields[1]))

            return interactions

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


class ReadGmt(ReadingData, ABC):
    """
    A class used to read the .gmt gene file inside pygna
    """

    def __init__(self, filename, read_descriptor=False):
        """
        :param filename: str, represents the path to the geneset file
        :param read_descriptor: bool, if the descriptor is given. Default = False
        """
        super().__init__(filename)
        self.filename = filename
        self.read_descriptor = read_descriptor

        self.gmt_data = self._ReadingData__readfile(self.filename)

    def _ReadingData__readfile(self, filename):
        """
        This method reads the geneset file into a variable
        :param filename: str, represents the path to the geneset file
        :returns gene_list: dict, represents the genes list
        """
        gene_lists = dict()
        with open(filename, "r") as f:
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


class ReadCsv(ReadingData, ABC):
    """
    A class used to read the .csv data file
    """

    def __init__(self, filename, sep=",", use_cols=None):
        """
        :param filename: str, represents the path to the data file
        :param sep: str, the separator to be used
        :param use_cols: list, columns used to be read and grouped
        """
        super().__init__(filename)
        self.filename = filename
        self.sep = sep
        self.use_cols = use_cols

        self.data = self._ReadingData__readfile()

    def _ReadingData__readfile(self):
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
