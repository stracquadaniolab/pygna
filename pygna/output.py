import logging
import networkx as nx
import os
import numpy as np
import statsmodels.stats.multitest as multi
import pandas as pd
import scipy
import tempfile
import shutil


class Output:
    """
    This class is used to print different data on files
    """

    def __init__(self, network_filename: str, output_table_results_file: str, analysis: str, geneset_file: str,
                 setnames: list, geneset_file_B: str = None, setnames_B: list = None):
        """
        :param network_filename: the file containing the network
        :param output_table_results_file: the output table that contains the results to use
        :param analysis: the type of analysis performed
        :param geneset_file: the geneset file use
        :param setnames: the names of the first geneset
        :param geneset_file_B: the second geneset file to use
        :param setnames_B: the names of the second geneset
        """
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

    def set_diffusion_matrix(self, diffusion_matrix_file: str) -> None:
        """
        Set the diffusion matrix file

        :param diffusion_matrix_file: set the diffusion matrix file to use

        Example
        ______
        >>> setnames = ["A", "B", "C"]
        >>> out = Output("networkfile.tsv", "results.csv", "myanalysis", "genset_a.csv", setnames)
        >>> out.set_diffusion_matrix("diffusion_matrix.csv")
        """
        self.diffusion_matrix_file = diffusion_matrix_file

    # Tables for stats
    def create_st_table_empirical(self) -> None:
        """
        Create the headings of the table that are going to be wrinnte in the csv file

        Example
        _______
        >>> setnames = ["A", "B", "C"]
        >>> out = Output("networkfile.tsv", "results.csv", "myanalysis", "genset_a.csv", setnames)
        >>> out.create_st_table_empirical()
        """
        tmp = tempfile.NamedTemporaryFile(mode='w+t', delete=False)
        self.table_file_name = tmp.name
        try:
            tmp.write("analysis,setname,n_mapped,n_geneset,number_of_permutations,observed,empirical_pvalue,mean(null),"
                      "var(null),network,geneset\n")
        finally:
            tmp.close()

    def close_temporary_table(self) -> None:
        """
        Remove the temporary file

        Example
        _______
        >>> setnames = ["A", "B", "C"]
        >>> out = Output("networkfile.tsv", "results.csv", "myanalysis", "genset_a.csv", setnames)
        >>> out.create_st_table_empirical()
        >>> out.close_temporary_table()
        """
        shutil.copy(self.table_file_name, self.output_table_results)
        os.remove(self.table_file_name)

    def update_st_table_empirical(self, setname: str, n_mapped: int, n_geneset: int, number_of_permutations: int,
                                  observed: int, empirical_pvalue: float, mean_null: np.mean, var_null: np.var) -> None:
        """
        Update the table content, adding a new line to the file

        :param setname: the name of the geneset
        :param n_mapped: the number of mapped genes
        :param n_geneset: the number of genesets
        :param number_of_permutations: the number of permutations
        :param observed: value of observed genes
        :param empirical_pvalue: value of the empirical p-value
        :param mean_null: mean of the null distribution
        :param var_null: var of the null distribution

        Example
        _______
        >>> setnames = ["A", "B", "C"]
        >>> out = Output("networkfile.tsv", "results.csv", "myanalysis", "genset_a.csv", setnames)
        >>> out.create_st_table_empirical()
        >>> out.update_st_table_empirical(setname, n_mapped, n_geneset, number_of_permutations, observed, pvalue=0.001, mean_null=np.mean(0.11), var_null=np.var(0.2))
        """
        setname = setname.replace(",", "_")
        with open(self.table_file_name, "a") as f:
            f.write(",".join([str(x) for x in
                              [self.analysis, setname, n_mapped, n_geneset, number_of_permutations, observed,
                               empirical_pvalue, mean_null, var_null, self.network_filename,
                               self.geneset_filename]]) + "\n")

    # Tables for comparisons
    def create_comparison_table_empirical(self) -> None:
        """
        Write the hadings for the comparison table

        Example
        _______
        >>> setnames = ["A", "B", "C"]
        >>> out = Output("networkfile.tsv", "results.csv", "myanalysis", "genset_a.csv", setnames)
        >>> out.create_comparison_table_empirical()
        """
        tmp = tempfile.NamedTemporaryFile(mode='w+t', delete=False)
        self.table_file_name = tmp.name
        try:
            tmp.write("analysis,setname_A,setname_B,n_geneset_A,n_mapped_A,n_geneset_B,n_mapped_B,n_overlaps,"
                      "number_of_permutations,observed,empirical_pvalue,mean(null),var(null),network\n ")
        finally:
            tmp.close()

    def update_comparison_table_empirical(self, setname_A: str, setname_B: str, n_geneset_A: int, n_mapped_A: int,
                                          n_geneset_B: int, n_mapped_B: int, n_overlaps: int,
                                          number_of_permutations: int, observed: int, empirical_pvalue: float,
                                          mean_null: np.mean, var_null: np.var) -> None:
        """
        Update the content of the comparison table, adding a new row on the file

        :param setname_A: the name of the geneset A
        :param setname_B: the name of the geneset B
        :param n_geneset_A: the number of genes in the geneset A
        :param n_mapped_A: the number of mapped genes in geneset A
        :param n_geneset_B: the number of genes in the geneset B
        :param n_mapped_B: the number of mapped genes in geneset B
        :param n_overlaps: the number of overlaps
        :param number_of_permutations: number of performed permutations
        :param observed: number of observed genes
        :param empirical_pvalue: value of the empirical pvalue
        :param mean_null: mean of the null distribution
        :param var_null: variance of the null distribution

        Example
        _______
        >>> import itertools
        >>> import pygna.command as cm
        >>> import pygna.reading_class as rc
        >>> import pygna.statistical_comparison as sc
        >>> geneset_a = rc.ReadGmt("genset_file").get_geneset("brca")
        >>> setnames = [key for key in geneset_a.keys()]
        >>> network = rc.ReadTsv("network_file.tsv").get_network()
        >>> distance_matrix_filename = "distance_matrix.tsv"
        >>> in_memory = True
        >>> network = nx.Graph(network.subgraph(max(nx.connected_components(network), key=len)))
        >>> sp_diz = {"nodes": cm.read_distance_matrix(distance_matrix_filename, in_memory=in_memory)[0],
        ...           "matrix": cm.read_distance_matrix(distance_matrix_filename, in_memory=in_memory)[1]}
        >>> st_comparison = sc.StatisticalComparison(sc.comparison_shortest_path, network, diz=sp_diz, n_proc=2)
        >>> out = Output("networkfile.tsv", "results.csv", "myanalysis", "genset_a.csv", setnames)
        >>> for pair in itertools.combinations(setnames, 2):
        ...     observed, pvalue, null_d, a_mapped, b_mapped = st_comparison.comparison_empirical_pvalue(set(geneset_a[pair[0]]), set(geneset_a[pair[1]]), max_iter=number_of_permutations)
        ...     out.update_comparison_table_empirical(pair[0], pair[1], len(set(geneset_a[pair[0]])), a_mapped, len(set(geneset_a[pair[1]])), b_mapped, n_overlaps, number_of_permutations, observed, pvalue, np.mean(null_d), np.var(null_d))
        """
        setname_A = setname_A.replace(",", "_")
        setname_B = setname_B.replace(",", "_")
        with open(self.table_file_name, "a") as f:
            f.write(",".join([str(x) for x in
                              [self.analysis, setname_A, setname_B, n_geneset_A, n_mapped_A, n_geneset_B, n_mapped_B,
                               n_overlaps, number_of_permutations, observed, empirical_pvalue, mean_null, var_null,
                               self.network_filename]]) + "\n")

    def add_GMT_entry(self, key: str, descriptor: str, gene_list: str) -> None:
        """
        Add a gmt entry in the GMT file

        :param key: the key name to store
        :param descriptor: the descriptor of the gene list
        :param gene_list: the gene list to write

        Example
        _______
        >>> geneset = rc.ReadGmt("geneset_file.csv").get_geneset("brca")
        >>> setnames = [key for key in geneset.keys()]
        >>> import pygna.reading_class as rc
        >>> network = rc.ReadTsv("network_file.tsv").get_network()
        >>> out = Output("networkfile.tsv", "results.csv", "myanalysis", "genset_a.csv", setnames)
        >>> for setname, item in geneset.items():
        ...     item = set(item)
        ...     module = nx.subgraph(network, item)
        ...     lcc = sorted(list(nx.connected_components(module)), key=len, reverse=True)[0]
        ...     out.add_GMT_entry("brca", "topology_module", lcc)
        """
        try:
            self.GMT_dict[key]
        except KeyError:
            self.GMT_dict[key] = {}
            self.GMT_dict[key]["descriptor"] = descriptor
            self.GMT_dict[key]["genes"] = gene_list

        else:
            logging.warning("Key Already Exists: " + str(key))

    def create_GMT_output(self, output_gmt: str) -> None:
        """
        Write the GMT line on the GMT file

        :param output_gmt: the GMT to print

        Example
        _______
        >>> geneset = rc.ReadGmt("geneset_file.csv").get_geneset("brca")
        >>> setnames = [key for key in geneset.keys()]
        >>> import pygna.reading_class as rc
        >>> network = rc.ReadTsv("network_file.tsv").get_network()
        >>> out = Output("networkfile.tsv", "results.csv", "myanalysis", "genset_a.csv", setnames)
        >>> for setname, item in geneset.items():
        ...     item = set(item)
        ...     module = nx.subgraph(network, item)
        ...     lcc = sorted(list(nx.connected_components(module)), key=len, reverse=True)[0]
        ...     out.add_GMT_entry("brca", "topology_module", lcc)
        >>> out.create_GMT_output("output_lcc.gmt")
        """
        self.output_gmt = output_gmt
        print_GMT(self.GMT_dict, self.output_gmt)


def print_GMT(gmt_dictionary: dict, output_file: str) -> None:
    """
    Save the dictionary on a GMT file

    :param gmt_dictionary: the dictionary containing the data
    :param output_file: the file to save the data

    Example
    _______
    >>> gmt_dict = {"key": "dict_sets"}
    >>> print_GMT(gmt_dict, "mygmt.gmt")
    """
    with open(output_file, "w") as f:
        f.write("")

    for key, dict_set in gmt_dictionary.items():
        with open(output_file, "a") as f:
            genes_dict = '\t'.join(map(str, dict_set["genes"]))
            f.write(str(key) + "\t" + str(dict_set["descriptor"]) + "\t" + genes_dict + "\n")


def apply_multiple_testing_correction(table_file: str, pval_col: str = "empirical_pvalue", method: str = "fdr_bh",
                                      threshold: float = 0.1) -> None:
    """
    Apply the multiple testing correction and save the file in a csv file

    :param table_file: the name of the file to read
    :param pval_col: the name column containing the empirical pvalue
    :param method: the correction method to use
    :param threshold: the threshold to use in the method

    Example
    _______
    >>> table_filename = "pygna_comparison_results.csv"
    >>> apply_multiple_testing_correction(table_filename, pval_col="empirical_pvalue", method="fdr_bh", threshold=0.1)
    """
    with open(table_file, "r+") as f:
        table = pd.read_csv(f)

    rejects, pval, k, bonf = multi.multipletests(table[pval_col].values, alpha=float(threshold), method=method)
    table["rejects"] = rejects
    table["bh_pvalue"] = pval
    table["k"] = k
    table["bonf"] = bonf

    table = table.sort_values(by="bh_pvalue")

    table.to_csv(table_file, index=False)


def write_graph_summary(graph: nx.Graph, output_file: str, net_name: str = None) -> None:
    """
    This function takes a graph as input and writes the network properties in a text file

    :param graph: the graph to print
    :param output_file: the name of the file to print
    :param net_name: the name of the network

    Example
    _______
    >>> import pygna.reading_class as rc
    >>> text_output = "My summary stats"
    >>> network = rc.ReadTsv("mynetwork.tsv").get_network()
    >>> write_graph_summary(network, text_output, "mynetwork.tsv")
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
