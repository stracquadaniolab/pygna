import logging
import pickle
import networkx as nx
import os
import sys
from datetime import datetime
import glob
import numpy as np
import statsmodels.stats.multitest as multi
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
import tempfile
import shutil


class Output:
    def __init__(
        self,
        network_filename,
        output_table_results_file,
        analysis,
        geneset_file,
        setnames,
        geneset_file_B=None,
        setnames_B=None,
    ):
        self.network_filename = network_filename
        self.analysis = analysis
        self.output_table_results = output_table_results_file
        self.output_gmt = None
        self.text = []
        self.geneset_filename = geneset_file
        self.setnames = setnames
        self.geneset_filename_B = geneset_file_B
        self.setnames_B = setnames_B
        self.diffusion_matrix_file = None
        self.GMT_dict = {}

        if not self.output_table_results.endswith('.csv'):
            logging.warning('The output table is saved as csv file, the name does not match the file extension')


    def set_diffusion_matrix(self, diffusion_matrix_file):
        self.diffusion_matrix_file = diffusion_matrix_file

    ## Tables for stats
    def create_st_table_empirical(self):

        tmp = tempfile.NamedTemporaryFile(mode='w+t', delete=False)
        self.table_file_name = tmp.name
        try:
            tmp.write(
                "analysis,setname,n_mapped,n_geneset,number_of_permutations,observed,empirical_pvalue,mean(null),var(null),network,geneset\n"
            )
        finally:
            tmp.close()
    
    def close_temporary_table(self):
        shutil.copy(self.table_file_name, self.output_table_results)   
        os.remove(self.table_file_name)


    def update_st_table_empirical(
        self,
        setname,
        n_mapped,
        n_geneset,
        number_of_permutations,
        observed,
        empirical_pvalue,
        mean_null,
        var_null,
    ):

        with open(self.table_file_name, "a") as f:
            f.write(
                ",".join(
                    [
                        str(x)
                        for x in [
                            self.analysis,
                            setname,
                            n_mapped,
                            n_geneset,
                            number_of_permutations,
                            observed,
                            empirical_pvalue,
                            mean_null,
                            var_null,
                            self.network_filename,
                            self.geneset_filename,
                        ]
                    ]
                )
                + "\n"
            )

    ## Tables for comparisons
    def create_comparison_table_empirical(self):
        tmp = tempfile.NamedTemporaryFile(mode='w+t', delete=False)
        self.table_file_name = tmp.name
        try:
            tmp.write(
                "analysis,setname_A,setname_B,n_geneset_A,n_mapped_A,n_geneset_B,n_mapped_B,n_overlaps,number_of_permutations,observed,empirical_pvalue,mean(null),var(null),network\n"
            )
        finally:
            tmp.close()

    def update_comparison_table_empirical(
        self,
        setname_A,
        setname_B,
        n_geneset_A,
        n_mapped_A,
        n_geneset_B,
        n_mapped_B,
        n_overlaps,
        number_of_permutations,
        observed,
        empirical_pvalue,
        mean_null,
        var_null,
    ):
        with open(self.table_file_name, "a") as f:
            f.write(
                ",".join(
                    [
                        str(x)
                        for x in [
                            self.analysis,
                            setname_A,
                            setname_B,
                            n_geneset_A,
                            n_mapped_A,
                            n_geneset_B,
                            n_mapped_B,
                            n_overlaps,
                            number_of_permutations,
                            observed,
                            empirical_pvalue,
                            mean_null,
                            var_null,
                            self.network_filename,
                        ]
                    ]
                )
                + "\n"
            )

    def add_GMT_entry(self, key, descriptor, gene_list):

        try:
            self.GMT_dict[key]
        except KeyError:
            self.GMT_dict[key] = {}
            self.GMT_dict[key]["descriptor"] = descriptor
            self.GMT_dict[key]["genes"] = gene_list
            
        else:
            logging.warning("Key Already Exists: " + str(key))

    def create_GMT_output(self, output_gmt):
        self.output_gmt = output_gmt
        print_GMT(self.GMT_dict, self.output_gmt)


def print_GMT(GMT_dictionary, output_file):

    with open(output_file, "w") as f:
        f.write("")

    for key, dict_set in GMT_dictionary.items():
        with open(output_file, "a") as f:
            f.write(
                str(key)
                + "\t"
                + str(dict_set["descriptor"])
                + "\t"
                + "\t".join(dict_set["genes"])
                + "\n"
            )


def apply_multiple_testing_correction(
    table_file, pval_col="empirical_pvalue", method="fdr_bh", threshold=0.1
        ):

    with open(table_file, "r+") as f:
        table = pd.read_csv(f)

    rejects, pval, k, bonf = multi.multipletests(
        table[pval_col].values, alpha=float(threshold), method=method
    )
    table["rejects"] = rejects
    table["bh_pvalue"] = pval
    table["k"] = k
    table["bonf"] = bonf

    table = table.sort_values(by="bh_pvalue")

    table.to_csv(table_file, index=False)


def write_graph_summary(graph, output_file, net_name=None):

    """
    This function takes a graph as input and writes the network
    properties in a text file
    """

    if not net_name:
        net_name = 'network'

    D = dict(nx.degree(graph))
    degree = np.array(list(dict(nx.degree(graph)).values()))

    n_nodes = nx.number_of_nodes(graph)
    n_edges = nx.number_of_edges(graph)

    degrees = {k: v for k, v in D.items()}
    degrees = sorted(degrees.items(), key=lambda kv: kv[1])

    density = (2 * n_edges) / ((n_nodes) * (n_nodes - 1))

    with open(output_file, "w") as file1:
        file1.write("Network Summary for %s " % str(net_name))
        file1.write("\n---------------------------------------------------\n")
        file1.write("\nInfo: " + nx.info(graph))
        file1.write("\nOther Properties::\n ")
        file1.write("\n\t- Density: " + str(density))
        file1.write("\n\t- min degree = " + str(np.min(degree)))
        file1.write("\n\t- max degree = " + str(np.max(degree)))
        file1.write("\n\t- median degree = " + str(np.median(degree)))
        file1.write("\n\t- degree mode = " + str(scipy.stats.mode(degree)))
        file1.write("\n\t- disconnected nodes = " + str(np.sum(degree == 0)))
        file1.write("\n\t- average clustering" + str(nx.average_clustering(graph)))


    largest_cc = nx.Graph(graph.subgraph(max(nx.connected_components(graph), key=len)))

    D = dict(nx.degree(largest_cc))
    degree = np.array(list(dict(nx.degree(largest_cc)).values()))

    n_nodes = nx.number_of_nodes(largest_cc)
    n_edges = nx.number_of_edges(largest_cc)

    degrees = {k: v for k, v in D.items()}
    degrees = sorted(degrees.items(), key=lambda kv: kv[1])

    density = (2 * n_edges) / ((n_nodes) * (n_nodes - 1))

    with open(output_file, "a") as file1:
        file1.write("\n\nLargest Connected Component Summary for %s \n " % str(net_name))
        file1.write(
            "-----------------------------------------------------------------\n"
        )
        file1.write("\nInfo: " + nx.info(largest_cc))
        file1.write("\nLargest Connected Component::\n")

        file1.write("\n\t- Density: " + str(density))
        file1.write("\n\t- min degree = " + str(np.min(degree)))
        file1.write("\n\t- max degree = " + str(np.max(degree)))
        file1.write("\n\t- median degree = " + str(np.median(degree)))
        file1.write("\n\t- degree mode = " + str(scipy.stats.mode(degree)))
        file1.write("\n\t- disconnected nodes = " + str(np.sum(degree == 0)))
        file1.write("\n\t- average clustering" + str(nx.average_clustering(largest_cc)))

