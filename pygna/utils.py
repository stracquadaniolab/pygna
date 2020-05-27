import yaml
import pygna.converters as pc
import pygna.elaborators as pe


class YamlConfig:
    def __init__(self):
        pass

    def write_config(self, data, filename):
        with open(filename, "w") as outfile:
            yaml.dump(data, outfile, default_flow_style=True)

    def load_config(self, filename):
        with open(filename, "r") as stream:
            try:
                config = yaml.load(stream)
                return config

            except yaml.YAMLError as exc:
                return print(exc)


def convert_gmt(gmt_file: "gmt file to be converted",
                output_gmt_file: "output file",
                conversion: "e2s or s2e",
                converter_map_filename: "tsv table used to convert gene names",
                entrez_col: "name of the entrez column" = "NCBI Gene ID",
                symbol_col: "name of the symbol column" = "Approved symbol"):
    pc.GmtToGmtEnriched(gmt_file, output_gmt_file, conversion, entrez_col, symbol_col, converter_map_filename)


def geneset_from_table(input_file: "input csv file",
                       setname: "name of the set",
                       output_gmt: "output gmt name" = None,
                       output_csv: "output csv name" = None,
                       name_column: "column with the names" = "Unnamed: 0",
                       filter_column: "column with the values to be filtered" = "padj",
                       alternative: "alternative to use for the filter, with less the filter is applied <threshold, "
                                    "otherwise >= threshold" = "less",
                       threshold: "threshold for the filter" = 0.01,
                       descriptor: "descriptor for the gmt file" = None):
    pc.CsvToGmt(input_file, setname, filter_column, alternative, threshold, output_gmt, output_csv, name_column,
                descriptor)


def filter_table(table: "input csv file",
                 filter_column: "column with the values to be filtered" = "padj",
                 alternative: "alternative to use for the filter, with less the filter is applied <threshold, "
                              "otherwise >= threshold" = "less",
                 threshold: "threshold for the filter" = 0.01):
    return pe.TableElaboration.filter_table(table, filter_column, alternative, threshold)


def generate_group_gmt(input_table: "table to get the geneset from",
                       output_gmt: "output gmt file",
                       name_col='Gene',
                       group_col='Cancer',
                       descriptor='cancer_genes'):
    pc.GroupGmt(input_table, output_gmt, name_col, group_col, descriptor)


def convert_csv(csv_file: "csv file where to add a name column",
                conversion: "e2s or s2e",
                original_name_col: "column name to be converted",
                new_name_col: "name of the new column with the converted names",
                geneset: "the geneset to convert",
                converter_map_filename: "tsv table used to convert gene names" = "entrez_name.tsv",
                output_file: "if none, table is saved in the same input file" = None,
                entrez_col: "name of the entrez column" = "NCBI Gene ID",
                symbol_col: "name of the symbol column" = "Approved symbol"):
    pc.CsvToCsvEnriched(csv_file, conversion, original_name_col, new_name_col, geneset, entrez_col, symbol_col,
                        converter_map_filename, output_file)
