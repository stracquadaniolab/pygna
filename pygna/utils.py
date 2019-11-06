"""Docstring1
"""

import logging
import pickle
import numpy as np
import networkx as nx

import pygna.statistical_test as st
import pygna.statistical_comparison as sc
import pygna.diagnostic as diagnostic
import pygna.painter as paint
import pandas as pd
import scipy
import time
import scipy.linalg.interpolative
from copy import copy, deepcopy
import itertools
import tables
import seaborn as sns
from matplotlib import pyplot as plt


import yaml
import pandas as pd
import pygna.parser as parser
import pygna.output as output
import logging


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


class Converter:
    def __init__(
        self,
        mapping_table_file="../../../primary_data/entrez_name.tsv",
        entrez_col="NCBI Gene ID",
        symbol_col="Approved symbol",
    ):
        with open(mapping_table_file, "r") as f:
            self.map_table = pd.read_table(f)

        self.map_table = self.map_table.fillna("0")
        self.map_table["Approved symbol"] = self.map_table[
            "Approved symbol"
        ].str.upper()
        self.map_table["Synonyms"] = self.map_table["Synonyms"].str.upper()
        self.map_table["Previous symbols"] = self.map_table[
            "Previous symbols"
        ].str.upper()
        self.entrez_column = entrez_col
        self.symbol_column = symbol_col

    def entrez2symbol(self, geneset):
        unknown_counter = 0
        geneset_symbol = []
        for i in geneset:
            name = self.map_table[self.map_table[self.entrez_column] == int(i)][
                self.symbol_column
            ].values.tolist()
            if len(name) > 0:
                geneset_symbol.append(str(name[0]))
            else:
                unknown_counter += 1
                geneset_symbol.append("<" + i + ">")
        if unknown_counter > 0:
            logging.warning(
                "%d/%d terms that couldn't be mapped" % (unknown_counter, len(geneset))
            )
        return geneset_symbol

    def symbol2entrez(self, geneset):
        geneset_entrez = []
        unknown_counter = 0
        for i in geneset:
            if type(i) != str:
                print(i)
                i = str(i)
            i = i.upper()
            name = self.map_table[self.map_table[self.symbol_column].str.upper() == i][
                self.entrez_column
            ].values.tolist()
            if len(name) > 0:
                geneset_entrez.append(str(int(name[0])))
            else:
                unknown_counter += 1
                geneset_entrez.append("<" + i + ">")

        if unknown_counter > 0:
            logging.warning(
                "%d/%d terms that couldn't be mapped" % (unknown_counter, len(geneset))
            )

        return geneset_entrez


def convert_gmt(
    gmt_file: "gmt file to be converted",
    output_gmt_file: "output file",
    conversion: "e2s or s2e",
    converter_map_filename: "tsv table used to convert gene names" = "../../../primary_data/entrez_name.tsv",
    entrez_col: "name of the entrez column" = "NCBI Gene ID",
    symbol_col: "name of the symbol column" = "Approved symbol",
    ):

    ''' name conversion table '''

    GMTparser = parser.GMTParser()
    genesets_dict = GMTparser.read(gmt_file, read_descriptor=True)

    converter = Converter(converter_map_filename, entrez_col, symbol_col)

    if conversion == "e2s":
        for key, dict in genesets_dict.items():
            genesets_dict[key]["genes"] = converter.entrez2symbol(dict["genes"])

    elif conversion == "s2e":
        for key, dict in genesets_dict.items():
            genesets_dict[key]["genes"] = converter.symbol2entrez(dict["genes"])
    else:
        logging.error("conversion type not understood")

    output.print_GMT(genesets_dict, output_gmt_file)


def geneset_from_table(
    input_file: "input csv file",
    setname: "name of the set",
    output_gmt: "output gmt name" = None,
    output_csv: "output csv name" = None,
    name_column: "column with the names" = "Unnamed: 0",
    filter_column: "column with the values to be filtered" = "padj",
    alternative: "alternative to use for the filter, with less the filter is applied <threshold, otherwise >= threshold" = "less",
    threshold: "threshold for the filter" = 0.01,
    descriptor: "descriptor for the gmt file" = None,
    ):

    """
    This function converts a csv file to a gmt allowing to filter the elements
    using the values of one of the columns. The user can specify the column used to
    retrieve the name of the objects and the filter condition.
    The output can be either a gmt with the names of the genes that pass the filter
    or a csv with the whole filtered table, otherwise both can be created.
    """

    if input_file.endswith(".csv"):
        with open(input_file, "r") as f:
            table = pd.read_csv(f)
    else:
        logging.error("only csv files supported")

    table[name_column]=table[name_column].fillna(0).apply(int).apply(str)
    threshold = float(threshold)

    table=clean_table(table, stat_col=filter_column)

    table = filter_table(
        table, filter_column=filter_column, alternative=alternative, threshold=threshold
    )

    if output_gmt:

        if descriptor == None:
            descriptor = input_file.split("/")[-1]

        gmt_dict = {}
        gmt_dict[setname] = {}
        gmt_dict[setname]["descriptor"] = descriptor
        gmt_dict[setname]["genes"] = []

        geneset = table.loc[:, name_column].values.tolist()

        logging.info("geneset=" + str(geneset))
        gmt_dict[setname]["genes"] = geneset

        if output_gmt.endswith(".gmt"):
            output.print_GMT(gmt_dict, output_gmt)
        else:
            logging.error("specify gmt output")

    if output_csv:

        if output_csv.endswith(".csv"):
            table.to_csv(output_csv, sep=",", index=False)
        else:
            logging.error("specify csv output")


def filter_table(
    table: "input csv file",
    filter_column: "column with the values to be filtered" = "padj",
    alternative: "alternative to use for the filter, with less the filter is applied <threshold, otherwise >= threshold" = "less",
    threshold: "threshold for the filter" = 0.01,
    ):

    """
    This function filters a table according to a filter rule.
    """

    threshold = float(threshold)

    try:
        if alternative == "less":
            table = table[table[filter_column] < threshold]
        else:
            table = table[table[filter_column] >= threshold]
    except:
        logging.error("error in filtering")

    return table


def convert_csv_names(
    csv_file: "csv file where to add a name column",
    conversion: "e2s or s2e",
    original_name_col: "column name to be converted",
    new_name_col: "name of the new column with the converted names",
    output_file: "if none, table is saved in the same input file" = None,
    converter_map_filename: "tsv table used to convert gene names" = "../../../primary_data/entrez_name.tsv",
    entrez_col: "name of the entrez column" = "NCBI Gene ID",
    symbol_col: "name of the symbol column" = "Approved symbol",
    ):

    with open(csv_file, "r") as f:
        table = pd.read_csv(f)

    converter = Converter(converter_map_filename, entrez_col, symbol_col)

    if conversion == "e2s":
        table[new_name_col] = converter.entrez2symbol(
            table[original_name_col].values.tolist()
        )

    elif conversion == "s2e":
        table[new_name_col] = converter.symbol2entrez(
            table[original_name_col].values.tolist()
        )

    else:
        logging.error("conversion type not understood")

    if not output_file:
        output_file = csv_file

    table.to_csv(output_file, index=False)


def clean_table(table, stat_col="stat"):
    logging.info("original table has %d rows" % len(table))
    table = table.dropna(subset=[stat_col])
    logging.info("cleaned table has %d rows" % len(table))
    return table

def generate_group_gmt( input_table:"table to get the geneset from",
                        output_gmt:'output_gmt_file',
                        name_col = 'Gene',
                        group_col = 'Cancer',
                        descriptor = 'cancer_genes',
                            ):

    '''
    This function generates a gmt file of multiple setnames.
    From the table file, it groups the names in the group_col (the column you
    want to use to group them) and prints the genes in the name_col.
    Set the descriptor according to your needs
    '''

    if input_table.endswith('.csv'):
        with open(input_table, 'r') as f:
            table = pd.read_csv(f, usecols = [name_col, group_col])
    elif (input_table.endswith('.tsv') or input_table.endswith('.txt')):        
        with open(input_table, 'r') as f:
            table = pd.read_csv(f, sep='\t', usecols = [name_col, group_col])
    else:
        sys.exit('pass correct input (csv/tsv/txt)')

    diz = {}
    for g, group in table.groupby([group_col]):
        if len(group)<10:
            print('warning: %s has less than 10 genes' %g)
        diz[g]={}
        diz[g]['genes']=group[name_col].astype(str).values.tolist()
        diz[g]['descriptor'] = descriptor
    
    out.print_GMT(diz, output_gmt)
