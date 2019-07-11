"""PyGNA commands main module
"""

from scipy.cluster.hierarchy import dendrogram, linkage, fcluster

import logging
import pickle
import numpy as np
from pygna.utils import YamlConfig
import pygna.output as out
import networkx as nx
import pygna.parser as ps
import pygna.statistical_test as st
import pygna.statistical_comparison as sc
import pygna.KS_test as KS
import pygna.diagnostic as diagnostic
import pygna.painter as paint
import pygna.utils as utils
import pygna.statistical_diffusion as sd
from pygna.utils import Converter
import pandas as pd
import scipy
import time
import scipy.linalg.interpolative
from copy import copy, deepcopy
import itertools
import tables
import seaborn as sns
from matplotlib import pyplot as plt


def __read_distance_matrix(distance_matrix_filename, in_memory=False):

    """ 
    Reads the large matrix .Uses a hdf5 file to work with it
    :params in_memory: pass the flag to store matrix in memory
    """

    if distance_matrix_filename.endswith("hdf5"):

        hdf5_nodes, hdf5_data = ps.Hdf5MatrixParser().read(
            distance_matrix_filename, in_memory
        )
        if type(hdf5_nodes[0]) == bytes:
            hdf5_nodes = [i.decode() for i in hdf5_nodes]
    else:
        logging.error("Invalid input format for matrix, use .hdf5")

    return hdf5_nodes, hdf5_data


################################################################################
######### SUMMARY ##############################################################
################################################################################


def network_summary(
    network_file, output_folder, output_name, geneset_file=None, setname=None
):
    """
    This function saves the principal info on a graph:
    - network properties
    - degree distribution

    If a geneset/setname is passed to the function, the properties of
    the subgraph are evaluated.
    """

    network = ps.__load_network(network_file)

    if geneset_file:
        geneset = ps.__load_geneset(geneset_file, setname)
        for setname, item in geneset.items():
            graph = nx.subgraph(network, item)
            out.write_graph_summary(graph, output_folder, output_name)

    else:

        out.write_graph_summary(network, output_folder, output_name)


################################################################################
######### SINGLE SET ANALYSES ##################################################
################################################################################


def test_topology_total_degree(
    network_file: "network file",
    geneset_file: "GMT geneset file",
    outpath: "output folder where to place all generated files",
    prefix: "prefix to add to all generated files",
    setname: "Geneset to analyse" = None,
    size_cut: "removes all genesets with a mapped length < size_cut" = 20,
    number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
    cores: "Number of cores for the multiprocessing" = 1,
    show_results: "barplot of results" = False,
    show_null: "plot null distribution" = False,
    symbol: "True if we want to print th output in symbol names" = False,
    ):

    """
        Performs the analysis of total degree of the .

        It computes a p-value for the ratio of total degree
        of the geneset being bigger than the one expected by chance
        for a geneset of the same size.
    """

    network = ps.__load_network(network_file)

    geneset = ps.__load_geneset(geneset_file, setname)

    setnames = [key for key in geneset.keys()]

    output1 = out.Output(
        network_file, outpath, prefix, "topology_total_degree", geneset_file, setnames
    )
    logging.info("Output-folder= " + output1.output)
    output1.create_st_table_empirical("table_total_degree")

    st_test = st.StatisticalTest(st.geneset_total_degree_statistic, network)

    for setname, item in geneset.items():

        if len(item) > size_cut:

            item = set(item)

            observed, pvalue, null_d, n_mapped, n_geneset = st_test.empirical_pvalue(
                item,
                max_iter=number_of_permutations,
                alternative="greater",
                cores=cores,
            )

            logging.info("Setname:" + setname)
            if n_mapped < size_cut:
                logging.info(
                    "%s remove from results since nodes mapped are < %d"
                    % (setname, size_cut)
                )
            else:
                logging.info("Observed: %g p-value: %g" % (observed, pvalue))
                logging.info(
                    "Null mean: %g null variance: %g"
                    % (np.mean(null_d), np.var(null_d))
                )

                output1.update_st_table_empirical(
                    setname,
                    n_mapped,
                    n_geneset,
                    number_of_permutations,
                    observed,
                    pvalue,
                    np.mean(null_d),
                    np.var(null_d),
                )
                if show_null:
                    diagnostic.plot_null_distribution(
                        null_d, observed, output1.output, setname=setname
                    )

    output1.save_output_summary()

    if show_results:
        paint.paint_datasets_stats(output1.output_table, output1.output, "total_degree")

def test_topology_internal_degree(
    network_file: "network file",
    geneset_file: "GMT geneset file",
    outpath: "output folder where to place all generated files",
    prefix: "prefix to add to all generated files",
    setname: "Geneset to analyse" = None,
    size_cut: "removes all genesets with a mapped length < size_cut" = 20,
    number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
    cores: "Number of cores for the multiprocessing" = 1,
    show_results: "barplot of results" = False,
    show_null: "plot null distribution" = False,
    symbol: "True if we want to print th output in symbol names" = False,
    ):

    """
        Performs the analysis of internal degree.

        It computes a p-value for the ratio of internal degree
        of the geneset being bigger than the one expected by chance
        for a geneset of the same size.
    """

    network = ps.__load_network(network_file)

    geneset = ps.__load_geneset(geneset_file, setname)

    setnames = [key for key in geneset.keys()]

    output1 = out.Output(
        network_file, outpath, prefix, "topology_internal_degree", geneset_file, setnames
    )
    logging.info("Output-folder= " + output1.output)
    output1.create_st_table_empirical("table_internal_degree")

    st_test = st.StatisticalTest(st.geneset_internal_degree_statistic, network)

    for setname, item in geneset.items():

        item = set(item)
        if len(item) > size_cut:
            observed, pvalue, null_d, n_mapped, n_geneset = st_test.empirical_pvalue(
                item,
                max_iter=number_of_permutations,
                alternative="greater",
                cores=cores,
            )

            logging.info("Setname:" + setname)

            if n_mapped < size_cut:
                logging.info(
                    "%s remove from results since nodes mapped are < %d"
                    % (setname, size_cut)
                )
            else:
                logging.info("Observed: %g p-value: %g" % (observed, pvalue))
                logging.info(
                    "Null mean: %g null variance: %g"
                    % (np.mean(null_d), np.var(null_d))
                )
                output1.update_st_table_empirical(
                    setname,
                    n_mapped,
                    n_geneset,
                    number_of_permutations,
                    observed,
                    pvalue,
                    np.mean(null_d),
                    np.var(null_d),
                )
                if show_null:
                    diagnostic.plot_null_distribution(
                        null_d, observed, output1.output, setname=setname
                    )

    output1.save_output_summary()

    if show_results:
        paint.paint_datasets_stats(
            output1.output_table, output1.output, "internal_degree"
        )

def test_topology_rwr(
    network_file: "network file, use a network with weights",
    geneset_file: "GMT geneset file",
    RW_dict_file: "hdf5 RWR matrix obtained with pygna ",
    outpath: "output folder where to place all generated files",
    prefix: "prefix to add to all generated files",
    setname: "Geneset to analyse" = None,
    size_cut: "removes all genesets with a mapped length < size_cut" = 20,
    weight: "RW" = "RW",
    number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
    cores: "Number of cores for the multiprocessing" = 1,
    in_memory: "set if you want the large matrix to be read in memory" = False,
    show_matrix: "plotting flag, if true the diffusion matrix for each geneset is saved " = False,
    show_results: "barplot of results" = False,
    show_null: "plot null distribution" = False,
    ):
    """
        Performs the analysis of random walk probabilities. 
        Given the RW matrix ( either normal random walk or RW with restart),
        it compares the probability of walking between the genes in the geneset 
        compared to those of walking between the nodes
        of a geneset with the same size
    """

    network = ps.__load_network(network_file)
    network = nx.Graph(network.subgraph(max(nx.connected_components(network), key=len)))

    geneset = ps.__load_geneset(geneset_file, setname)
    RW_dict = {}  

    RW_dict["nodes"], RW_dict["matrix"] = __read_distance_matrix(
        RW_dict_file, in_memory=in_memory
    )

    setnames = [key for key in geneset.keys()]
    output1 = out.Output(
        network_file, outpath, prefix, "topology_rwr", geneset_file, setnames
    )

    logging.info("Output-folder= " + output1.output)
    output1.create_st_table_empirical("table_RW")

    output1.set_diffusion_matrix(RW_dict_file)

    st_test = st.StatisticalTest(st.geneset_RW_statistic, network, RW_dict)

    for setname, item in geneset.items():

        item = set(item).intersection(set(list(network.nodes)))
        if len(item) > size_cut:
            # test
            observed, pvalue, null_d, n_mapped, n_geneset = st_test.empirical_pvalue(
                item, max_iter=number_of_permutations, alternative="greater"
            )
            
            if len(item) > 0 and show_null:
                logging.info("Plotting diagnostic" + str(output1.output))
                diagnostic.plot_null_distribution(
                    null_d, observed, output1.output, setname=setname
                )

            logging.info("Setname:" + setname)
            logging.info("Observed: %g p-value: %g" % (observed, pvalue))
            logging.info(
                "Null mean: %g null variance: %g" % (np.mean(null_d), np.var(null_d))
            )
            # saving output
            output1.update_st_table_empirical(
                setname,
                n_mapped,
                n_geneset,
                number_of_permutations,
                observed,
                pvalue,
                np.mean(null_d),
                np.var(null_d),
            )

            logging.info("show matrix: " + str(show_matrix))
            if len(item) > 0 and show_matrix:
                print(len(item))
                print(item)
                painter_rw = paint.Painter_RW(
                    network, output1.output, setname, RW_dict, item
                )
                logging.info("Painting matrix")
                painter_rw.plot_matrix()

        else:
            logging.info(
                "%s removed from results since nodes mapped are < %d"
                % (setname, size_cut)
            )

    output1.save_output_summary()

    if show_results:
        paint.paint_datasets_stats(output1.output_table, output1.output, "RW")

def test_topology_module(
    network_file: "network file",
    geneset_file: "GMT geneset file",
    outpath: "output folder where to place all generated files",
    prefix: "prefix to add to all generated files",
    setname: "Geneset to analyse" = None,
    size_cut: "removes all genesets with a mapped length < size_cut" = 20,
    number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
    cores: "Number of cores for the multiprocessing" = 1,
    create_output_LCC: "flag for creating a GMT file with the LCC lists" = False,
    show_results: "barplot of results" = False,
    show_null: "plot null distribution" = False,
    symbol: "True if we want to print th output in symbol names" = False,
    ):
    """
        Performs geneset network topology module analysis.

        It computes a p-value for the largest connected component
        of the geneset being bigger than the one expected by chance
        for a geneset of the same size.
    """
    network = ps.__load_network(network_file)

    geneset = ps.__load_geneset(geneset_file, setname)

    setnames = [key for key in geneset.keys()]

    output1 = out.Output(
        network_file, outpath, prefix, "topology_module", geneset_file, setnames
    )
    logging.info("Output-folder= " + output1.output)
    output1.create_st_table_empirical("table_module")

    st_test = st.StatisticalTest(st.geneset_module_statistic, network)

    for setname, item in geneset.items():

        item = set(item)
        if len(item) > size_cut:
            if create_output_LCC:

                module = nx.subgraph(network, item)
                if len(module.nodes) > 0:
                    LCC = sorted(
                        list(nx.connected_components(module)), key=len, reverse=True
                    )[0]
                else:
                    LCC = []
                if symbol:
                    converter = Converter()
                    LCC = converter.entrez2symbol(LCC)

                output1.add_GMT_entry(setname, "topology_module", LCC)

            observed, pvalue, null_d, n_mapped, n_geneset = st_test.empirical_pvalue(
                item, max_iter=number_of_permutations, alternative="greater"
            )
            # observed_z=(observed-np.mean(null_d))/np.std(null_d)

            logging.info("Setname:" + setname)
            if n_mapped < size_cut:
                logging.info(
                    "%s remove from results since nodes mapped are < %d"
                    % (setname, size_cut)
                )
            else:
                logging.info("Observed: %g p-value: %g" % (observed, pvalue))
                logging.info(
                    "Null mean: %g null variance: %g"
                    % (np.mean(null_d), np.var(null_d))
                )

                output1.update_st_table_empirical(
                    setname,
                    n_mapped,
                    n_geneset,
                    number_of_permutations,
                    observed,
                    pvalue,
                    np.mean(null_d),
                    np.var(null_d),
                )
                if show_null:
                    diagnostic.plot_null_distribution(
                        null_d, observed, output1.output, setname=setname
                    )

    output1.save_output_summary()

    if create_output_LCC:
        output1.create_GMT_output()

    if show_results:
        paint.paint_datasets_stats(output1.output_table, output1.output, "module")

def test_topology_sp(
    network_file: "network file",
    geneset_file: "GMT geneset file",
    distance_matrix_filename: "distance hdf5 matrix file generated by pygna",
    outpath: "output folder where to place all generated files",
    prefix: "prefix to add to all generated files",
    setname: "Geneset to analyse" = None,
    size_cut: "removes all genesets with a mapped length < size_cut" = 20,
    number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
    cores: "Number of cores for the multiprocessing" = 1,
    in_memory: "set if you want the large matrix to be read in memory" = False,
    show_matrix: "plotting flag, if true the distance matrix for each geneset is saved " = False,
    show_results: "barplot of results" = False,
    ):
    """
        Performs geneset network topology shortest path analysis.

        It computes a p-value for the average shortest path length
        of the geneset being smaller than expected by chance
        for a geneset of the same size.

    """

    network = ps.__load_network(network_file)
    network = nx.Graph(network.subgraph(max(nx.connected_components(network), key=len)))

    geneset = ps.__load_geneset(geneset_file, setname)

    diz = {}
    diz["nodes"], diz["matrix"] = __read_distance_matrix(
        distance_matrix_filename, in_memory=in_memory
    )

    setnames = [key for key in geneset.keys()]

    output1 = out.Output(
        distance_matrix_filename,
        outpath,
        prefix,
        "topology_sp",
        geneset_file,
        setnames,
    )
    logging.info("Output-folder= " + output1.output)
    output1.create_st_table_empirical("table_SP")

    st_test = st.StatisticalTest(st.geneset_localisation_statistic, network, diz)

    for setname, item in geneset.items():

        item = set(item).intersection(set(list(network.nodes)))
        if len(item) > size_cut:
            observed, pvalue, null_d, n_mapped, n_geneset = st_test.empirical_pvalue(
                item, max_iter=number_of_permutations
            )

            logging.info("Observed: %g p-value: %g" % (observed, pvalue))
            logging.info(
                "Null mean: %g null variance: %g" % (np.mean(null_d), np.var(null_d))
            )
            logging.info("1th percentile: %g " % (np.percentile(null_d, 1)))

            output1.update_st_table_empirical(
                setname,
                n_mapped,
                n_geneset,
                number_of_permutations,
                observed,
                pvalue,
                np.mean(null_d),
                np.var(null_d),
            )
        else:
            logging.info(
                "%s remove from results since nodes mapped are < %d"
                % (setname, size_cut)
            )

    output1.save_output_summary()

    if show_results:
        paint.paint_datasets_stats(output1.output_table, output1.output, "SP")


################################################################################
######### Degree Distribution test #############################################
################################################################################


def test_degree_distribution(
    network_file: "network file",
    geneset_file: "GMT geneset file",
    outpath: "output folder where to place all generated files",
    prefix: "prefix to add to all generated files",
    setname: "Geneset to analyse" = None,
    size_cut: "removes all genesets with a mapped length < size_cut" = 20,
    show_results: "barplot of results" = False,
    ):
    """
        Performs degree distribution test.
        Kolmogorov-Smirnov statistic on 2 samples.
        H0 is that the geneset is drawn from the same distribution of all the other nodes. 
        H0 rejected if statistic is greater.
    """
    network = ps.__load_network(network_file)

    geneset = ps.__load_geneset(geneset_file, setname)

    setnames = [key for key in geneset.keys()]

    output1 = out.Output(
        network_file, outpath, prefix, "test_degree", geneset_file, setnames
    )
    logging.info("Output-folder= " + output1.output)
    output1.create_st_table_empirical("table_degree")

    st_test = KS.KSTest(KS.degree_distribution, network)

    for setname, item in geneset.items():

        item = set(item)
        if len(item) > size_cut:
            observed, pvalue, n_mapped, n_geneset = st_test.apply_test(item)

            logging.info("Setname:" + setname)
            if n_mapped < size_cut:
                logging.info(
                    "%s remove from results since nodes mapped are < %d"
                    % (setname, size_cut)
                )
            else:
                logging.info("Observed: %g p-value: %g" % (observed, pvalue))

                output1.update_st_table_empirical(
                    setname, n_mapped, n_geneset, 1, observed, pvalue, 0, 0
                )

    output1.save_output_summary()

    if show_results:
        paint.paint_datasets_stats(output1.output_table, output1.output, "test_degree")


################################################################################
######### Diffusion test #######################################################
################################################################################


def test_diffusion_weights(
    network_file: "network file, use a network with weights",
    geneset_file: "csv geneset file",
    RW_dict_file: "pickle obtained calculating a random walk matrix",
    outpath: "output folder where to place all generated files",
    prefix: "prefix to add to all generated files",
    name_column: "Column to use as name" = "stat",
    weight_column: "Column to use as weight (default is deseq)" = "stat",
    filter_column: "Column used to define the significant genes" = "padj",
    filter_condition: "Condition for significance" = "less",
    filter_threshold: "threshold for significance" = 0.01,
    weight: "RW" = "RW",
    number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
    cores: "Number of cores for the multiprocessing" = 1,
    show_matrix: "plotting flag, if true the diffusion matrix for each geneset is saved " = False,
    show_results: "barplot of results" = False,
    show_null: "plot null distribution" = False,
    ):

    """
        Performs the analysis of random walk diffusion weights. Given the RW matrix
        ( either normal random walk or RW with restart), it applies the weights from a differential
        expression file (or a general csv). Diffusion weights for the geneset are used as statistic
        and compared to an empirical null distribution generated by permutations
    """

    network = ps.__load_network(network_file)
    network = nx.Graph(network.subgraph(max(nx.connected_components(network), key=len)))

    with open(geneset_file, "r") as f:
        table = pd.read_csv(geneset_file, sep=",")
    if len(table.columns) < 2:
        logging.warning(
            "Error: the function takes a csv file as input, the read file has less than 2 columns, check that the table is comma separated"
        )

    table = utils.clean_table(table, stat_col=weight_column)
    geneset = utils.filter_table(
        table,
        filter_column=filter_column,
        alternative=filter_condition,
        threshold=filter_threshold,
    )[name_column]
    print(geneset)

    RW_dict = {}
    RW_dict["nodes"], RW_dict["matrix"] = __read_distance_matrix(RW_dict_file)

    # Decoding nodes
    nodes = [i.decode() for i in RW_dict["nodes"]]

    # setting output
    output1 = out.Output(
        network_file,
        outpath,
        prefix,
        "test_diffusion_weights",
        geneset_file,
        geneset_file,
    )
    logging.info("Output-folder= " + output1.output)
    output1.create_st_table_empirical("table_diffusion_weights")
    output1.set_diffusion_matrix(RW_dict_file)

    # initialising test

    st_test = sd.DiffusionTest(
        sd.weights_diffusion_statistic,
        nodes,
        RW_dict["matrix"],
        table,
        names_col=name_column,
        weights_col=weight_column,
    )

    observed, pvalue, null_d, n_mapped, n_geneset = st_test.empirical_pvalue(
        geneset, max_iter=number_of_permutations, alternative="greater", cores=cores
    )
    # observed_z=(observed-np.mean(null_d))/np.std(null_d

    if show_null:
        logging.info("Plotting diagnostic %s" % (str(output1.output)))

        diagnostic.plot_null_distribution(
            null_d, observed, output1.output, setname="diffusion"
        )

    logging.info("Observed: %g p-value: %g" % (observed, pvalue))
    logging.info("Null mean: %g null variance: %g" % (np.mean(null_d), np.var(null_d)))
    output1.update_st_table_empirical(
        geneset_file,
        n_mapped,
        n_geneset,
        number_of_permutations,
        observed,
        pvalue,
        np.mean(null_d),
        np.var(null_d),
    )

    output1.save_output_summary()


################################################################################
######### associations and COMPARISONS #########################################
################################################################################

def test_association_sp(
    network_file: "network file",
    A_geneset_file: "GMT geneset file, if it's the only parameter passed the analysis is gonna be run on all the couples of datasets, otherwise specify the other files and setnames",
    distance_matrix_filename: "distance matrix file generated by pygna",
    outpath: "output folder where to place all generated files",
    prefix: "prefix to add to all generated files",
    setname_A: "Geneset A to analyse" = None,
    B_geneset_file: "GMT geneset file" = None,
    setname_B: "Geneset B to analyse" = None,
    size_cut: "removes all genesets with a mapped length < size_cut" = 20,
    keep: "if true, keeps the geneset B unpermuted" = False,
    cores: "Number of cores for the multiprocessing" = 1,
    in_memory: "set if you want the large matrix to be read in memory" = False,
    number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
    show_matrix: "plotting flag, if true the diffusion matrix for each pair of genesets is saved" = False,
    show_results: "heatmap of results" = False,
    ):
    """
        Performs comparison of network location analysis. If the flag
        --keep  is passed, the B geneset is kept
        fixed, and doesnt't get permuted

        It computes a p-value for the shortest path distance
        between two genesets being smaller than expected by chance

        If only A_geneset_file is passed the analysis is run on all the couples
        of sets in the file,
        if both A_geneset_file and B_geneset_file are passed, one can specify
        the setnames for both, if there is only one geneset in the file, setname_X
        can be omitted,
        if both sets are in the same file, B_geneset_file can be not specified,
        but setnames are needed.
    """

    if keep:
        analysis_name_str = "association_sp"
    else:
        analysis_name_str = "comparison_sp"

    network = ps.__load_network(network_file)
    network = nx.Graph(network.subgraph(max(nx.connected_components(network), key=len)))

    sp_diz = {}
    sp_diz["nodes"], sp_diz["matrix"] = __read_distance_matrix(
        distance_matrix_filename, in_memory=in_memory
    )

    if setname_A and setname_B == None and B_geneset_file == None:
        logging.error(" this analysis requires at least two genesets ")

    geneset_A = ps.__load_geneset(A_geneset_file, setname_A)

    if B_geneset_file:
        geneset_B = ps.__load_geneset(B_geneset_file, setname_B)
    else:
        if setname_B:
            geneset_B = ps.__load_geneset(A_geneset_file, setname_B)
        else:
            geneset_B = None

    st_comparison = sc.StatisticalComparison(
        sc.comparison_shortest_path, network, diz=sp_diz, n_proc=cores
    )

    if not geneset_B:  # Analysis of genesets inside a single file

        logging.info("Analysising all the sets in " + A_geneset_file)
        setnames = [key for key in geneset_A.keys()]
        logging.info("Setnames: " + str(setnames))

        # Creating the output table
        output1 = out.Output(
            network_file, outpath, prefix, analysis_name_str, A_geneset_file, setnames
        )
        output1.add_output_text(" distance matrix = " + str(distance_matrix_filename))
        logging.info("Output-folder= " + output1.output)
        output1.create_comparison_table_empirical("table_" + analysis_name_str)

        for pair in itertools.combinations(setnames, 2):
            logging.info("Analysing " + str(pair[0]) + " and " + str(pair[1]))
            overlaps = set(geneset_A[pair[0]]).intersection(set(geneset_A[pair[1]]))
            logging.info("There are %d genes shared between A and B" % len(overlaps))

            observed, pvalue, null_d, A_mapped, B_mapped = st_comparison.comparison_empirical_pvalue(
                set(geneset_A[pair[0]]),
                set(geneset_A[pair[1]]),
                max_iter=number_of_permutations,
                keep=keep,
            )
            # Save the results
            output1.update_comparison_table_empirical(
                pair[0],
                pair[1],
                len(set(geneset_A[pair[0]])),
                A_mapped,
                len(set(geneset_A[pair[1]])),
                B_mapped,
                len(overlaps),
                number_of_permutations,
                observed,
                pvalue,
                np.mean(null_d),
                np.var(null_d),
            )

    else:  # Analysis of genesets into two different gmt files

        logging.info("geneset_A contains %d sets", (len(geneset_A)))
        sets_A = [key for key in geneset_A.keys()]
        logging.info("Setnames in A: " + str(sets_A))
        logging.info("geneset_B contains %d sets", (len(geneset_B)))
        sets_B = [key for key in geneset_B.keys()]
        logging.info("Setnames in B: " + str(sets_B))

        output1 = out.Output(
            network_file,
            outpath,
            prefix,
            analysis_name_str,
            A_geneset_file,
            sets_A,
            B_geneset_file,
            sets_B,
        )
        output1.add_output_text(" distance matrix = " + str(distance_matrix_filename))
        logging.info("Output-folder= " + output1.output)
        output1.create_comparison_table_empirical("table_" + analysis_name_str)

        for set_A, item_A in geneset_A.items():
            for set_B, item_B in geneset_B.items():
                n_overlaps = len(set(item_A).intersection(set(item_B)))

                observed, pvalue, null_d, A_mapped, B_mapped = st_comparison.comparison_empirical_pvalue(
                    set(item_A), set(item_B), max_iter=number_of_permutations, keep=keep
                )

                logging.info("Observed: %g p-value: %g" % (observed, pvalue))
                logging.info(
                    "Null mean: %g null variance: %g"
                    % (np.mean(null_d), np.var(null_d))
                )
                logging.info("1th percentile: %g " % (np.percentile(null_d, 1)))

                output1.update_comparison_table_empirical(
                    set_A,
                    set_B,
                    len(set(item_A)),
                    A_mapped,
                    len(set(item_B)),
                    B_mapped,
                    n_overlaps,
                    number_of_permutations,
                    observed,
                    pvalue,
                    np.mean(null_d),
                    np.var(null_d),
                )

    output1.save_output_summary()  # Save Summary of Analysis

    if show_results:
        paint.paint_comparison_stats(output1.output_table, output1.output + "/", "sp")


def test_association_rwr(
    network_file: "network file",
    A_geneset_file: "GMT geneset file",
    RW_dict_file: ".hdf5 file with the RWR matrix obtained by pygna",
    outpath: "output folder where to place all generated files",
    prefix: "prefix to add to all generated files",
    setname_A: "Geneset A to analyse" = None,
    B_geneset_file: "GMT geneset file" = None,
    setname_B: "Geneset B to analyse" = None,
    size_cut: "removes all genesets with a mapped length < size_cut" = 20,
    keep: "if true, keeps the geneset B unpermuted" =False,
    cores: "Number of cores for the multiprocessing" = 1,
    in_memory: "set if you want the large matrix to be read in memory" = False,
    number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
    show_matrix: "plotting flag, if true the diffusion matrix for each pair of genesets is saved " = False,
    show_results: "heatmap of results" = False,
):
    """
        Performs comparison of network location analysis.

        It computes a p-value for the shortest path distance
        between two genesets being smaller than expected by chance

        If only A_geneset_file is passed the analysis is run on all the couples
        of sets in the file,
        if both A_geneset_file and B_geneset_file are passed, one can specify
        the setnames for both, if there is only one geneset in the file, setname_X
        can be omitted,
        if both sets are in the same file, B_geneset_file can be not specified,
        but setnames are needed.
    """

    if keep:
        analysis_name_str = "association_rwr"
    else:
        analysis_name_str = "comparison_rwr"

    network = ps.__load_network(network_file)
    network = nx.Graph(network.subgraph(max(nx.connected_components(network), key=len)))

    # Read datasets
    if setname_A and setname_B == None and B_geneset_file == None:
        logging.error(" this analysis requires at least two genesets ")

    logging.info("Reading diffusion matrix")
    RW_dict = {}
    RW_dict["nodes"], RW_dict["matrix"] = __read_distance_matrix(
        RW_dict_file, in_memory=in_memory
    )

    geneset_A = ps.__load_geneset(A_geneset_file, setname_A)

    if B_geneset_file:
        geneset_B = ps.__load_geneset(B_geneset_file, setname_B)
    else:
        if setname_B:
            geneset_B = ps.__load_geneset(A_geneset_file, setname_B)
        else:
            geneset_B = None

    st_comparison = sc.StatisticalComparison(
        sc.comparison_random_walk, network, n_proc=cores, diz=RW_dict
    )

    if not geneset_B:

        logging.info("Analysising all the sets in " + A_geneset_file)
        setnames = [key for key in geneset_A.keys()]
        logging.info("Setnames: " + str(setnames))

        output1 = out.Output(
            network_file, outpath, prefix, analysis_name_str, A_geneset_file, setnames
        )
        output1.add_output_text(" RW matrix = " + str(RW_dict_file))
        logging.info("Output-folder= " + output1.output)
        output1.create_comparison_table_empirical("table_" + analysis_name_str)

        for pair in itertools.combinations(setnames, 2):
            if len(pair[0]) > size_cut and len(pair[1]) > size_cut:

                logging.info("Analysing " + str(pair[0]) + " and " + str(pair[1]))

                n_overlaps = len(
                    set(geneset_A[pair[0]]).intersection(set(geneset_A[pair[1]]))
                )

                observed, pvalue, null_d, A_mapped, B_mapped = st_comparison.comparison_empirical_pvalue(
                    set(geneset_A[pair[0]]),
                    set(geneset_A[pair[1]]),
                    max_iter=number_of_permutations,
                    alternative="greater",
                    keep=keep,
                )

                output1.update_comparison_table_empirical(
                    pair[0],
                    pair[1],
                    len(set(geneset_A[pair[0]])),
                    A_mapped,
                    len(set(geneset_A[pair[1]])),
                    B_mapped,
                    n_overlaps,
                    number_of_permutations,
                    observed,
                    pvalue,
                    np.mean(null_d),
                    np.var(null_d),
                )

    else:

        logging.info("geneset_A contains %d sets" % ((len(geneset_A))))
        sets_A = [key for key in geneset_A.keys()]
        logging.info("Setnames in A: " + str(sets_A))
        logging.info("geneset_B contains %d sets" % ((len(geneset_B))))
        sets_B = [key for key in geneset_B.keys()]
        logging.info("Setnames in B: " + str(sets_B))

        output1 = out.Output(
            network_file,
            outpath,
            prefix,
            analysis_name_str,
            A_geneset_file,
            sets_A,
            B_geneset_file,
            sets_B,
        )
        logging.info("Output-folder= " + output1.output)
        output1.create_comparison_table_empirical("table_" + analysis_name_str)
        output1.add_output_text(" RW matrix = " + str(RW_dict_file))

        for set_A, item_A in geneset_A.items():
            for set_B, item_B in geneset_B.items():

                if len(item_A) > size_cut and len(item_B) > size_cut:
                    logging.info("Analysing " + str(set_A) + " and " + str(set_B))
                    n_overlaps = len(set(item_A).intersection(set(item_B)))

                    observed, pvalue, null_d, A_mapped, B_mapped = st_comparison.comparison_empirical_pvalue(
                        set(item_A),
                        set(item_B),
                        max_iter=number_of_permutations,
                        alternative="greater",
                        keep=keep,
                    )

                    logging.info("Observed: %g p-value: %g" % (observed, pvalue))
                    logging.info(
                        "Null mean: %g null variance: %g"
                        % (np.mean(null_d), np.var(null_d))
                    )
                    logging.info("1th percentile: %g " % (np.percentile(null_d, 1)))

                    output1.update_comparison_table_empirical(
                        set_A,
                        set_B,
                        len(set(item_A)),
                        A_mapped,
                        len(set(item_B)),
                        B_mapped,
                        n_overlaps,
                        number_of_permutations,
                        observed,
                        pvalue,
                        np.mean(null_d),
                        np.var(null_d),
                    )

    output1.save_output_summary()  # Save Summary of Analysis
    if show_results:
        sg= setname_B == None
        paint.paint_comparison_RW(output1.output_table, output1.output, "RWR", single_geneset=sg)


################################################################################
######### BUILDING FUNCTIONS ###################################################
################################################################################


def build_distance_matrix(
    network_file: "network file",
    output_file: "distance matrix output file, use a format between .lm.txt and .hdf5",
    giant_component_only: "compute the shortest paths only for nodes in the giant component" = True,
):
    """
        Build a shortest path distance matrix for a given network.
        Matrix can be saved as a lm.txt file or a .hdf5 one.
    """

    # network = ps.__load_network(network_file)
    network = ps.__load_network(network_file)
    if giant_component_only:
        network = network.subgraph(max(nx.connected_components(network), key=len))

    print(type(network))
    distance_matrix = nx.all_pairs_shortest_path_length(network)
    print(type(distance_matrix))

    nodes = list(network.nodes())

    if output_file.endswith(".hdf5"):
        with tables.open_file(output_file, mode="w") as hdf5_file:

            # create a hdf5 file with two objects:
            # - one is the nodes array,
            hdf5_nodes = hdf5_file.create_array(hdf5_file.root, "nodes", nodes)
            # -  the other is the shorthest path distance matrix
            hdf5_data = hdf5_file.create_array(
                hdf5_file.root, "matrix", np.zeros((len(nodes), len(nodes)))
            )

            for node_i, k in distance_matrix:
                i = nodes.index(node_i)

                for node_j, sp in k.items():
                    j = nodes.index(node_j)
                    if j >= i:  # saves only the upper triangular matrix
                        hdf5_data[i, j] = sp

    elif output_file.endswith(".lm.txt"):
        # creates the lm.txt file .
        # The format is node_i \t node_j \t shortest_path
        visited_nodes = []
        with open(output_file, "w") as f:
            f.write(str(nodes)[1:-1].replace(",", "\t") + "\n")
            for node_i, k in distance_matrix:
                i = nodes.index(node_i)
                print(i)
                visited_nodes.append(node_i)
                for node_j, sp in k.items():
                    j = nodes.index(node_j)
                    if node_j not in visited_nodes:
                        f.write(str(i) + "\t" + str(j) + "\t" + str(sp) + "\n")
    else:
        logging.error("enter a valid extension for the output file")

    # pickle.dump(distance_matrix, open(output_file, 'wb'))


def build_RWR_diffusion(
    network_file: "network file",
    beta=0.85,
    output_file: "distance matrix output file (use .hdf5) " = None,
    output_figure: "diffusion matrix output filename" = None,
):
    """
        Build the RWR_diffusion_matrix
    """

    network = ps.__load_network(network_file)
    network = network.subgraph(max(nx.connected_components(network), key=len))

    nodes = list(network.nodes())

    logging.info("Beginning to calculate RWR matrix")

    A = nx.adjacency_matrix(network)
    K = 1 / np.array(list(dict(network.degree()).values()))
    D = scipy.sparse.dia_matrix((K, [0]), shape=(A.shape))
    A = A.dot(D)
    n = np.shape(A)[1]

    if output_file.endswith(".hdf5"):
        with tables.open_file(output_file, mode="w") as hdf5_file:

            # create a hdf5 file with two objects:
            # - one is the nodes array,
            hdf5_nodes = hdf5_file.create_array(hdf5_file.root, "nodes", nodes)
            # -  the other is the RWR matrix
            hdf5_data = hdf5_file.create_array(
                hdf5_file.root,
                "matrix",
                beta * np.linalg.inv(np.eye(n) - (1.0 - beta) * A),
            )

            logging.info("Saving network")
    else:
        return beta * np.linalg.inv(np.eye(n) - (1.0 - beta) * A)

    if output_figure:
        RW_dict = {}
        RW_dict["nodes"], RW_dict["matrix"] = __read_distance_matrix(output_file)
        diagnostic.plot_diffusion_matrix(RW_dict["nodes"], RW_dict["matrix"], output_figure)


def build_graph(
    network_file: "network file",
    geneset_file: "geneset file",
    output_folder: "graphml network for visualisation",
    setname: "setname" = None,
    giant_component_only: "compute the shortest paths only for nodes in the giant component" = True,
):
    """
        Build a shortest path distance matrix for a given network.
    """
    network = ps.__load_network(network_file)
    geneset = ps.__load_geneset(geneset_file, setname)

    if giant_component_only:
        network = network.subgraph(max(nx.connected_components(network), key=len))

    for setname in geneset:
        print(setname)
        new_network = nx.Graph()
        new_network_minimal = nx.Graph()
        for s in geneset[setname]:
            if s in network.nodes():
                l0 = np.inf
                for t in geneset[setname]:
                    if (t in network.nodes()) & (s != t):
                        path = nx.shortest_path(network, source=s, target=t)
                        if l0 > len(path):
                            l0 = len(path)
                            new_network_minimal.add_path(path)
                        new_network.add_path(path)

        dict_nodes = {}
        for n in new_network.nodes():
            if n in geneset[setname]:
                dict_nodes[n] = True
            else:
                dict_nodes[n] = False
        nx.set_node_attributes(new_network, dict_nodes, "in_subset")
        nx.write_graphml(new_network, output_folder + setname + ".graphml")

        dict_nodes = {}
        for n in new_network_minimal.nodes():
            if n in geneset[setname]:
                dict_nodes[n] = True
            else:
                dict_nodes[n] = False
        nx.set_node_attributes(new_network_minimal, dict_nodes, "in_subset")

        nx.write_graphml(
            new_network_minimal, output_folder + setname + "_minimal.graphml"
        )


def network_graphml(
    network_file: "network file",
    geneset_file: "geneset file",
    output_folder: "graphml network for visualisation",
    prefix: "prefix for the new file",
    setname: "setname" = None,
    giant_component_only: "compute the shortest paths only for nodes in the giant component" = True,
):
    """
        Build a shortest path distance matrix for a given network.
    """
    network = ps.__load_network(network_file)
    geneset = ps.__load_geneset(geneset_file, setname)

    if giant_component_only:
        network = network.subgraph(max(nx.connected_components(network), key=len))

    for setname in geneset:
        # print(setname)
        # new_network = nx.Graph()
        # new_network_minimal = nx.Graph()
        # for s in geneset[setname]:
        #     if s in network.nodes():
        #         l0 = np.inf
        #         for t in geneset[setname]:
        #             if (t in network.nodes()) & (s != t):
        #                 path = nx.shortest_path(network, source=s, target=t)
        #                 if l0 > len(path):
        #                     l0 = len(path)
        #                     new_network_minimal.add_path(path)
        #                 new_network.add_path(path)

        dict_nodes = {}
        for n in network.nodes():
            if n in geneset[setname]:
                dict_nodes[n] = True
            else:
                dict_nodes[n] = False
        nx.set_node_attributes(network, dict_nodes, setname)

    nx.write_graphml(network, output_folder + prefix + ".graphml")
