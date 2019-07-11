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


class Output:
    def __init__(
        self,
        network_filename,
        output_path,
        prefix,
        analysis,
        geneset_file,
        setnames,
        geneset_file_B=None,
        setnames_B=None,
    ):

        self.network_filename = network_filename
        self.analysis = analysis
        self.output_path = output_path
        self.output = self.output_path + prefix + "_"
        self.text = []
        self.geneset_filename = geneset_file
        self.setnames = setnames
        self.geneset_filename_B = geneset_file_B
        self.setnames_B = setnames_B
        self.diffusion_matrix_file = None
        self.GMT_dict = {}

        try:
            # today = datetime.now()
            os.listdir(self.output_path)

            # dirs=glob.glob(self.output_path+today.strftime('%Y%m%d')+"_*")
            # print(dirs)
            # if dirs:
            #    iter=[int(i[len(self.output_path+today.strftime('%Y%m%d'))+1:]) for i in dirs]
            #    print(iter)
            #    self.output_folder=today.strftime('%Y%m%d')+"_"+str(max(iter)+1)
            # else:
            #    self.output_folder= today.strftime('%Y%m%d')+"_0"
            # os.mkdir(self.output_path+self.output_folder)

        except FileNotFoundError:
            logging.error("Output path doesn't exists")
            sys.exit(-1)

    def set_diffusion_matrix(self, diffusion_matrix_file):
        self.diffusion_matrix_file = diffusion_matrix_file

    def add_output_text(self, text):

        """Add text to the output. text an be both a string or a list
        of values convertible to strings.
        """

        if type(text) == str:
            self.text.append(text)
        elif type(text) == list:
            for t in text:
                self.text.append(str(t))
        else:
            logging.error("Text needs to be a string or a list of strings")

    def save_output_summary(self):

        """Summary.txt is written using the input configurations
        and the text that has been added to the output instance"""

        with open(self.output + "summary.txt", "w") as file1:
            file1.write("Network= " + str(self.network_filename))
            file1.write("\n Input file= " + str(self.geneset_filename))
            file1.write("\n Analysis= " + str(self.analysis))
            file1.write("\n Setnames = " + str(self.setnames))
            if self.geneset_filename_B:
                file1.write("\n Geneset file B= " + str(self.geneset_filename_B))
            if self.setnames_B:
                file1.write("\n Setname B= " + str(self.setnames_B))
            if self.diffusion_matrix_file:
                file1.write("\n Diffusion matrix= " + str(self.diffusion_matrix_file))
            for line in self.text:
                file1.write("\n" + line)

    ## Tables for stats
    def create_st_table_empirical(self, output_table_file):

        self.output_table = self.output + output_table_file + ".csv"
        with open(self.output_table, "w") as f:
            f.write(
                "analysis,setname,n_mapped,n_geneset,number_of_permutations,observed,empirical_pvalue,mean(null),var(null),network,geneset\n"
            )

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

        with open(self.output_table, "a") as f:
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
    def create_comparison_table_empirical(self, output_table_file):
        self.output_table = self.output + output_table_file + ".csv"
        with open(self.output_table, "w") as f:
            f.write(
                "analysis,setname_A,setname_B,n_geneset_A,n_mapped_A,n_geneset_B,n_mapped_B,n_overlaps,number_of_permutations,observed,empirical_pvalue,mean(null),var(null),network\n"
            )

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
        with open(self.output_table, "a") as f:
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
            logging.info("Key added to dictionary" + str(key))
        else:
            logging.info("Key Already Exists: " + str(key))

    def create_GMT_output(self):
        output_file = self.output + "LCC_gene_list.gmt"
        print_GMT(self.GMT_dict, output_file)


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


def write_graph_summary(graph, output_folder, prefix):

    """
    This function takes a graph as input and writes the network properties in a file
    """

    D = dict(nx.degree(graph))
    degree = np.array(list(dict(nx.degree(graph)).values()))

    n_nodes = nx.number_of_nodes(graph)
    n_edges = nx.number_of_edges(graph)

    degrees = {k: v for k, v in D.items()}
    degrees = sorted(degrees.items(), key=lambda kv: kv[1])

    density = (2 * n_edges) / ((n_nodes) * (n_nodes - 1))

    with open(output_folder + prefix + "_graph_summary.txt", "w") as file1:
        file1.write("Network Summary for %s " % str(prefix))
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

    fig, axes = plt.subplots(1, figsize=(10, 10))
    g1 = sns.distplot(degree, hist=True, kde=True, rug=False, ax=axes)
    perc=np.percentile(degree,99)
    g1 = sns.distplot(degree[degree>perc], hist=False, kde=False, rug=True,color='r', ax=axes)
    sns.despine(ax=axes, top=True, bottom=False, right=True, left=True)
    g1.set_ylabel("Density")
    g1.set_xlabel("Degree")
    fig.savefig(output_folder + prefix + "_degree.pdf", format="pdf")
    fig.savefig(output_folder + prefix + "_degree.png", format="png")

    largest_cc = nx.Graph(graph.subgraph(max(nx.connected_components(graph), key=len)))

    D = dict(nx.degree(largest_cc))
    degree = np.array(list(dict(nx.degree(largest_cc)).values()))

    n_nodes = nx.number_of_nodes(largest_cc)
    n_edges = nx.number_of_edges(largest_cc)

    degrees = {k: v for k, v in D.items()}
    degrees = sorted(degrees.items(), key=lambda kv: kv[1])

    density = (2 * n_edges) / ((n_nodes) * (n_nodes - 1))

    with open(output_folder + prefix + "_graph_summary.txt", "a") as file1:
        file1.write("\n\nLargest Connected Component Summary for %s \n " % str(prefix))
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

