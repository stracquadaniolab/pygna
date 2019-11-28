from pygna.utilities import Utils
import logging


class TableElaboration(Utils):
    """
    This class contains methods to use to elaborate a table (a pd.dataframe)
    """

    @staticmethod
    def clean_table(table, stat_col="stat"):
        """
        This method clean the table from the N/A values
        :param table: pd.dataframe, representing the table to be cleaned
        :param stat_col: str, the column to be cleaned
        :return: pd.dataframe, the table cleaned
        """
        logging.info("original table has %d rows" % len(table))
        table = table.dropna(subset=[stat_col])
        logging.info("cleaned table has %d rows" % len(table))
        return table

    @staticmethod
    def filter_table(table, filter_column="padj", alternative="less", threshold=0.01):
        """
        This method filters a table according to a filter rule
        :param table: pd.dataframe, the table to be filtered
        :param filter_column: str, column with the values to be filtered
        :param alternative: str, alternative to use for the filterK with "less" the filter is applied <threshold;
        otherwise >= threshold
        :param threshold: float, threshold for the filter
        :return: pd.dataframe, the table filtered
        """
        threshold = float(threshold)
        try:
            if alternative == "less":
                table = table[table[filter_column] < threshold]
            else:
                table = table[table[filter_column] >= threshold]
        except Exception as e:
            logging.error("type error: " + str(e))
        return table
