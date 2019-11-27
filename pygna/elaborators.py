from pygna.utilities import Utils


class TableElaboration(Utils):
    def clean_table(table, stat_col="stat"):
        logging.info("original table has %d rows" % len(table))
        table = table.dropna(subset=[stat_col])
        logging.info("cleaned table has %d rows" % len(table))
        return table
