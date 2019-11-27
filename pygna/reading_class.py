from abc import ABC, abstractmethod

import pandas as pd


class ReadingData(ABC):
    """Abstract class used to read different types of file. Each subclass must implement the 'readfile' method"""
    def __init__(self, filename):
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
        :param int_type: Unknown
        """
        super().__init__(filename)
        self._filename = filename
        self._int_type = int_type
        self._pd_table = pd_table

        if not self._pd_table:
            self._interactions = self._ReadingData__readfile(self._filename)

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
                if self._int_type:
                    types = fields[3].split(";")
                    if self._int_type in types:
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
        return self._interactions

    def get_pd_data(self):
        """
        Returns the data into a pandas object
        :return: pd.dataframe, represents the data into a pandas object
        """
        return pd.read_table(self._filename)


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
        self._filename = filename
        self._read_descriptor = read_descriptor

        self._gmt_data = self._ReadingData__readfile(self._filename)

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
                if self._read_descriptor:
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
        return self._gmt_data


class ReadCsv(ReadingData, ABC):
    """
    A class used to read the .csv data file
    """

    def __init__(self, filename, sep=","):
        """
        :param filename: str, represents the path to the data file
        """
        super().__init__(filename)
        self._filename = filename
        self._sep = sep

        self._data = self._ReadingData__readfile(self._filename, sep)

    def _ReadingData__readfile(self, filename, sep):
        """
        This method read the file and saves the data inside a class attribute
        :param filename:  str, represents the path to the file
        :return: pd.dataframe, represents teh data read inside the .csv
        """
        with open(filename, "r") as f:
            table = pd.read_csv(f, sep=sep)
            return table

    def get_data(self):
        """
        Returns the data of the csv file
        :return: pd.dataframe, represents teh data read inside the .csv
        """
        return self._data

    def fill_na_column(self, name_column="gene_name"):
        """
        Fill the N/A values with a (str) 0
        :param name_column: name of the column to filter
        """
        self._data[name_column] = self._data[name_column].fillna(0).apply(int).apply(str)
