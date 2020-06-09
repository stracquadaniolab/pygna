import logging
import pandas as pd


class TableElaboration:
    """
    This class contains methods to use to elaborate a table
    """

    @staticmethod
    def clean_table(table: pd.DataFrame(), stat_col: str = "stat") -> pd.DataFrame():
        """
        This function clean the table from the N/A values

        :param table: dataframerepresenting the table to be cleaned
        :param stat_col: the column to be cleaned
        :return: the table cleaned from the N/A values
        """
        logging.info("original table has %d rows" % len(table))
        table = table.dropna(subset=[stat_col])
        logging.info("cleaned table has %d rows" % len(table))
        return table

    @staticmethod
    def filter_table(table: pd.DataFrame(), filter_column: str = "padj", alternative: str = "less",
                     threshold: float = 0.01) -> pd.DataFrame():
        """
        This method filters a table according to a filter rule

        :param table: The table to be filtered
        :param filter_column: Column with the values to be filtered
        :param alternative: Alternative to use for the filterK with "less" the filter is applied <threshold; otherwise >= threshold
        :param threshold: Threshold for the filter
        :return: The table filtered
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
