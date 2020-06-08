from pygna.elaborators import TableElaboration as tE
import pygna.reading_class as rc
import logging
import sys
import pandas as pd

import pygna.output as output


class Converters:
    """
    This class is has utilities that can be used to convert the data
    """

    def __init__(self):
        super(Converters, self).__init__()

    @classmethod
    def convert_e2s(cls, geneset: pd.DataFrame, tsv_data: pd.DataFrame, entrez_col: str = "NCBI Gene ID",
                    symbol_col: str = "Approved symbol") -> list:
        """
        Method to convert the entrez 2 string

        :param tsv_data: the dataframe to work on
        :param symbol_col: the column containing the symbols
        :param entrez_col: the column containing the entrez ID
        :param geneset: column containing the entrez to convert
        :return: list containing the string names
        """
        logging.info("Converting Entrez ID -> Symbols")
        unknown_counter = 0
        geneset_symbol = []
        for i in geneset:
            name = tsv_data[tsv_data[entrez_col] == int(i)][symbol_col].values.tolist()
            if len(name) > 0:
                geneset_symbol.append(str(name[0]))
            else:
                unknown_counter += 1
                geneset_symbol.append("<" + i + ">")
        if unknown_counter > 0:
            logging.warning("%d/%d terms that couldn't be mapped" % (unknown_counter, len(geneset)))

        return geneset_symbol

    @classmethod
    def convert_s2e(cls, geneset: pd.DataFrame, tsv_data: pd.DataFrame, entrez_col: str = "NCBI Gene ID",
                    symbol_col: str = "Approved symbol") -> list:
        """
        Method to convert the string 2 entrez

        :param tsv_data: the dataframe to work on
        :param symbol_col: the column containing the symbols
        :param entrez_col: the column containing the entrez ID
        :param geneset: column containing the strings to convert
        :return: list containing the entrez names
        """
        logging.info("Converting Symbols -> Entrez ID")
        geneset_entrez = []
        unknown_counter = 0
        for i in geneset:
            if type(i) != str:
                i = str(i)
            i = i.upper()
            name = tsv_data[tsv_data[symbol_col].str.upper() == i][entrez_col].values.tolist()
            if len(name) > 0:
                geneset_entrez.append(str(int(name[0])))
            else:
                unknown_counter += 1
                geneset_entrez.append("<" + i + ">")
        if unknown_counter > 0:
            logging.warning("%d/%d terms that couldn't be mapped" % (unknown_counter, len(geneset)))

        return geneset_entrez

    @staticmethod
    def _gmt_output(gmt_data: dict, gmt_output_file: str) -> None:
        """
        Save the dictionary on a GMT file

        :param gmt_data: the dictionary containing the data
        :param gmt_output_file: the file to save the data
        """
        output.print_GMT(gmt_data, gmt_output_file)


class CsvToCsvEnriched(Converters):
    """
    Class that is used to add a column with the entrezID or Symbols to a CSV file
    """

    def __init__(self, csv_file: pd.DataFrame, conversion: str, original_name_col: str, new_name_col: str, geneset: str,
                 entrez_col: str, symbol_col: str, converter_map_filename: str = "entrez_name.tsv",
                 output_file: str = None):
        """
        :param csv_file: dataframe with the data
        :param conversion: could be "e2s"-> Entrez2Symbols or "s2e" -> Symbol2Entrez
        :param original_name_col:  the column where to find the information to convert
        :param new_name_col: the name of the column that is going to contain the information
        :param geneset: the geneset to convert
        :param converter_map_filename:  the path to the .tsv used to convert the genes name
        :param output_file: [optional] the path of the output file
        :param entrez_col: the name of the entrez column
        :param symbol_col:  the name of the symbol column
        """
        super().__init__()
        logging.info("Adding values to the CSV file")
        self.filename = csv_file
        self.conversion = conversion
        self.original_name_col = original_name_col
        self.new_name_col = new_name_col
        self.geneset = geneset
        self.converter_map_file = converter_map_filename
        self.output = output_file
        self.entrez_col = entrez_col
        self.symbol_col = symbol_col

        self.map_table = rc.ReadTsv(self.converter_map_file, pd_table=True).get_data()
        self.file_data = rc.ReadCsv(self.filename).get_data()
        self._clean_table()

        if self.conversion == "e2s":
            self.file_data[self.new_name_col] = \
                super().convert_e2s(self.file_data[self.original_name_col].values.tolist(),
                                    self.map_table, self.entrez_col, self.symbol_col)
        elif self.conversion == "s2e":
            self.file_data[self.new_name_col] = \
                super().convert_s2e(self.file_data[self.original_name_col].values.tolist(),
                                    self.map_table, self.entrez_col, self.symbol_col)
        else:
            logging.error("Conversion type not understood")

        if self.output:
            self._csv_output()

    def get_data(self) -> pd.DataFrame:
        """
        Return the conversion result

        :return: dataframe with the e2s or s2e added as column
        """
        return self.file_data

    def _csv_output(self) -> None:
        """
        Print the output to a csv
        """
        output_file = self.output
        self.filename.to_csv(output_file, index=False)

    def _clean_table(self) -> None:
        """
        Method to make all upper and clean the table from null values
        """
        self.map_table = self.map_table.fillna("0")
        self.map_table["Approved symbol"] = self.map_table["Approved symbol"].str.upper()
        self.map_table["Synonyms"] = self.map_table["Synonyms"].str.upper()
        self.map_table["Previous symbols"] = self.map_table["Previous symbols"].str.upper()


class CsvToGmt(Converters):
    """
    This Class converts a csv file to a gmt allowing to filter the elements using the values of one of the columns.
    The user can specify the column used to retrieve the name of the objects and the filter condition. The output
    can be either a gmt with the names of the genes that pass the filter or a csv with the whole filtered table,
    otherwise both can be created.
    """

    def __init__(self, input_file: str, setname: str, filter_column: str, alternative: str, threshold: float,
                 output_gmt: str = None, output_csv: str = None, name_column: str = "Unnamed: 0",
                 descriptor: str = None):
        """
        :param input_file: str, the csv file
        :param setname: str, the name of the set
        :param output_gmt: str, output gmt name
        :param output_csv: str, output csv name
        :param name_column: str, column with the names
        :param filter_column: str, column with the values to be filtered
        :param alternative: str, alternative to use for the filterK with "less" the filter is applied <threshold; otherwise >= threshold
        :param threshold: float, threshold for the filter
        :param descriptor: str, descriptor for the gmt file
        """
        super().__init__()
        self.input_file = input_file
        self.setname = setname
        self.output_gmt = output_gmt
        self.output_csv = output_csv
        self.name_column = name_column
        self.filter_column = filter_column
        self.alternative = alternative
        self.threshold = threshold
        self.descriptor = descriptor

        self.table = rc.ReadCsv(self.input_file, column_to_fill=self.name_column).get_data()
        self.table = self._elaborate()

        if self.output_gmt:
            self._process_gmt()

        if self.output_csv:
            self._csv_output()

    def _elaborate(self) -> pd.DataFrame:
        """
        This method performs the cleaning and the filtering of the table

        :return: dataframe representing the cleaned and filter table
        """
        table = tE.clean_table(self.table, self.filter_column)
        table = tE.filter_table(table, filter_column=self.filter_column, alternative=self.alternative,
                                threshold=self.threshold)
        return table

    def _process_gmt(self) -> None:
        """
        This method parse the results and save them in a GMT file
        """
        if self.descriptor is None:
            self.descriptor = self.input_file.split("/")[-1]
        gmt_dict = {self.setname: {}}
        gmt_dict[self.setname]["descriptor"] = self.descriptor
        gmt_dict[self.setname]["genes"] = []
        geneset = self.table.loc[:, self.name_column].values.tolist()

        logging.info("geneset=" + str(geneset))
        gmt_dict[self.setname]["genes"] = geneset

        if self.output_gmt.endswith(".gmt"):
            super()._gmt_output(gmt_dict, self.output_gmt)
        else:
            logging.error("specify gmt output")

    def _csv_output(self) -> None:
        """
        This method save the pd.dataframe in a CSV file
        """

        if self.output_csv.endswith(".csv"):
            self.table.to_csv(self.output_csv, sep=",", index=False)
        else:
            logging.error("specify csv output")


class GmtToGmtEnriched(Converters):
    """
    This Class converts a GMT file, adding information about the Entrez ID or the symbol
    """

    def __init__(self, gmt_file: str, output_gmt_file: str, conversion: str, entrez_col: str, symbol_col: str,
                 converter_map_filename: str = "entrez_name.tsv"):
        """
        :param gmt_file: the input GMT file path
        :param output_gmt_file: the output GMT file path
        :param conversion: could be "e2s"-> Entrez2Symbols or "s2e" -> Symbol2Entrez
        :param entrez_col: the name of the entrez column
        :param symbol_col: the name of the symbol column
        :param converter_map_filename: the path to the .tsv used to convert the genes name
        """
        super().__init__()
        self.gmt_file = gmt_file
        self.output_gmt_file = output_gmt_file
        self.conversion = conversion
        self.entrez_col = entrez_col
        self.symbol_col = symbol_col
        self.converter_map_filename = converter_map_filename

        self.gmt_data = rc.ReadGmt(self.gmt_file, True).get_data()
        self.tsv_data = rc.ReadTsv(self.converter_map_filename, pd_table=True).get_data()

        if self.conversion == "e2s":
            for k, d in self.gmt_data.items():
                self.gmt_data[k]["genes"] = super().convert_e2s(d["genes"], self.tsv_data,
                                                                self.entrez_col, self.symbol_col)
        elif self.conversion == "s2e":
            for k, d in self.gmt_data.items():
                self.gmt_data[k]["genes"] = super().convert_s2e(d["genes"], self.tsv_data,
                                                                self.entrez_col, self.symbol_col)
        else:
            logging.error("Conversion type not understood")
        super()._gmt_output(self.gmt_data, self.output_gmt_file)


class GroupGmt(Converters):
    """
    This function generates a gmt file of multiple setnames. From the table file, it groups the names in the
    group_col (the column you want to use to group them) and prints the genes in the name_col. Set the descriptor
    according to your needs
    """

    def __init__(self, input_table: str, output_gmt: str, name_col: str = "Gene", group_col: str = "Cancer",
                 descriptor: str = "cancer_genes"):
        """
        :param input_table: the filename path
        :param output_gmt: the output gmt file path
        :param name_col: the name of the column to write the genes
        :param group_col:the name of the column to group
        :param descriptor: the descriptor to use
        """
        super().__init__()
        self.input_table = input_table
        self.output_gmt = output_gmt
        self.name_col = name_col
        self.group_col = group_col
        self.descriptor = descriptor

        if self.input_table.endswith(".csv"):
            self.table = rc.ReadCsv(self.input_table, use_cols=[self.name_col, self.group_col]).get_data()
        elif self.input_table.endswith("tsv") or self.input_table.endswith("txt"):
            self.table = rc.ReadCsv(self.input_table, sep="\t", use_cols=[self.name_col, self.group_col]).get_data()
        else:
            sys.exit('Pass correct input (csv/tsv/txt)')
        self.gmt_data = self._elaborate()
        super()._gmt_output(self.gmt_data, self.output_gmt)

    def _elaborate(self) -> dict:
        """
        This method elaborates the table and returns a dictionary with the grouped columns

        :return: with the genes and the descriptor
        """
        diz = {}
        for g, group in self.table.groupby([self.group_col]):
            if len(group) < 10:
                print('warning: %s has less than 10 genes' % g)
            diz[g] = {}
            diz[g]['genes'] = group[self.name_col].astype(str).values.tolist()
            diz[g]['descriptor'] = self.descriptor
        return diz
