from pygna.utilities import Utils
import pygna.reading_class as rc
import logging


class Converters(Utils):
    """
    Class representing a Converter utility
    """

    def __init__(self):
        super(Converters, self).__init__()

    @classmethod
    def convert_e2s(cls, geneset, tsv_data, entrez_col="NCBI Gene ID", symbol_col="Approved symbol"):
        """
        Method to convert the entrez 2 string
        :param tsv_data: pd.dataframe, the dataframe to work on
        :param symbol_col: str, the column containing the symbols
        :param entrez_col: str, the column containing the entrez ID
        :param geneset: pd.dataframe, column containing the entrez to convert
        :return: list, containing the string names
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
    def convert_s2e(cls, geneset, tsv_data, entrez_col="NCBI Gene ID", symbol_col="Approved symbol"):
        """
        Method to convert the string 2 entrez
        :param tsv_data: pd.dataframe, the dataframe to work on
        :param symbol_col: str, the column containing the symbols
        :param entrez_col: str, the column containing the entrez ID
        :param geneset: pd.dataframe, column containing the strings to convert
        :return: list, containing the entrez names
        """
        logging.info("Converting Symbols -> Entrez ID...")
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


class ConvertCsvNames(Converters):
    """
    Class that is used to add a column with the entrezID or Symbols to a CSV file
    """

    def __init__(self, csv_file, conversion, original_name_col, new_name_col, geneset, entrez_col,
                 symbol_col, converter_map_filename="entrez_name.tsv", output_file=None):
        """
        :param csv_file: str, representing the path of the file
        :param conversion: str, could be "e2s"-> Entrez2Symbols or "s2e" -> Symbol2Entrez
        :param original_name_col: str, the column where to find the information to convert
        :param new_name_col: str, the name of the column that is going to contain the information
        :param geneset: str, the geneset to convert
        :param converter_map_filename: str, the path to the .tsv used to convert the genes name
        :param output_file: [optional] str, the path of the output file
        :param entrez_col: str, the name of the entrez column
        :param symbol_col: str, the name of the symbol column
        """
        super().__init__()
        logging.info("Adding values to the CSV file")
        self._filename = csv_file
        self._conversion = conversion
        self._original_name_col = original_name_col
        self._new_name_col = new_name_col
        self._geneset = geneset
        self._converter_map_file = converter_map_filename
        self._output = output_file
        self._entrez_col = entrez_col
        self._symbol_col = symbol_col

        self._map_table = rc.ReadTsv(self._converter_map_file, pd_table=True).get_pd_data()
        self._file_data = rc.ReadCsv(self._filename).get_data()
        self._clean_table()

        if self._conversion == "e2s":
            self._file_data[self._new_name_col] = \
                Converters.convert_e2s(self._file_data[self._original_name_col].values.tolist(),
                                       self._map_table, self._entrez_col, self._symbol_col)
        elif self._conversion == "s2e":
            self._file_data[self._new_name_col] = \
                Converters.convert_s2e(self._file_data[self._original_name_col].values.tolist(),
                                       self._map_table, self._entrez_col, self._symbol_col)
        else:
            logging.error("Conversion type not understood")

        if self._output:
            self._print_output()

    def get_data(self):
        """
        Return the conversion result
        :return: pd.dataframe object with the e2s or s2e added as column
        """
        return self._file_data

    def _print_output(self):
        """
        Method to print the output to file
        :return: null
        """
        output_file = self._output
        self._filename.to_csv(output_file, index=False)

    def _clean_table(self):
        """
        Method to make all upper and clean the table from null values
        :return: null
        """
        self._map_table = self._map_table.fillna("0")
        self._map_table["Approved symbol"] = self._map_table["Approved symbol"].str.upper()
        self._map_table["Synonyms"] = self._map_table["Synonyms"].str.upper()
        self._map_table["Previous symbols"] = self._map_table["Previous symbols"].str.upper()


class CsvToGmtFiltered(Converters):
    """
    This Class converts a csv file to a gmt allowing to filter the elements using the values of one of the columns.
    The user can specify the column used to retrieve the name of the objects and the filter condition. The output
    can be either a gmt with the names of the genes that pass the filter or a csv with the whole filtered table,
    otherwise both can be created.
    """

    def __init__(self, input_file, setname, output_gmt=None, output_csv=None, name_column="Unnamed: 0",
                 filter_column="padj", alternative="less", threshold=0.01, descriptor=None):
        """
        :param input_file: str, the csv file
        :param setname: str, the name of the set
        :param output_gmt: str, output gmt name
        :param output_csv: str, output csv name
        :param name_column: str, column with the names
        :param filter_column: str, column with the values to be filtered
        :param alternative: str, alternative to use for the filterK with "less" the filter is applied <threshold;
        otherwise >= threshold
        :param threshold: float, threshold for the filter
        :param descriptor: str, descriptor for the gmt file
        """
        super().__init__()
        self._input_file = input_file
        self._setname = setname
        self._output_gmt = output_gmt
        self._output_csv = output_csv
        self._name_column = name_column
        self._filter_column = filter_column
        self._alternative = alternative
        self._threshold = float(threshold)
        self._descriptor = descriptor

        self._table = rc.ReadCsv(self._input_file).fill_na_column(self._name_column).get_data()
        self._table = self._elaborate()

        if self._output_gmt:
            self._gmt_output()

        if self._output_csv:
            self._csv_output()

    def _elaborate(self):
        """
        This method performs the cleaning and the filtering of the table
        :return: pd.dataframe, representing the cleaned and filter table
        """
        table = clean_table(self._table, stat_col=self._filter_column)
        table = filter_table(
            table, filter_column=self._filter_column, alternative=self._alternative, threshold=self._threshold
        )
        return table

    def _gmt_output(self):
        """
        This method parse the results and save them on a GMT file
        :return: null
        """
        if self._descriptor is None:
            self._descriptor = self._input_file.split("/")[-1]
        gmt_dict = {self._setname: {}}
        gmt_dict[self._setname]["descriptor"] = self._descriptor
        gmt_dict[self._setname]["genes"] = []
        geneset = self._table.loc[:, self._name_column].values.tolist()

        logging.info("geneset=" + str(geneset))
        gmt_dict[self._setname]["genes"] = geneset

        # TODO Maybe it's better to automatically add an extension
        if self._output_gmt.endswith(".gmt"):
            # TODO Fix below
            output.print_GMT(gmt_dict, self._output_gmt)
        else:
            logging.error("specify gmt output")

    def _csv_output(self):
        """
        This method sabe the pd.dataframe on a CSV file
        :return: null
        """
        # TODO Maybe it's better to automatically add an extension
        if self._output_csv.endswith(".csv"):
            self._table.to_csv(self._output_csv, sep=",", index=False)
        else:
            logging.error("specify csv output")
