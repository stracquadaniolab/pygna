"""Docstring1
"""

from scipy.cluster.hierarchy import dendrogram, linkage, fcluster

import logging
import pickle
import numpy as np
from geneset_network_analysis_tool.utils import YamlConfig
import geneset_network_analysis_tool.output as out
import networkx as nx
import geneset_network_analysis_tool.parser as ps
import geneset_network_analysis_tool.statistical_test as st
import geneset_network_analysis_tool.statistical_comparison as sc
import geneset_network_analysis_tool.KS_test as KS
import geneset_network_analysis_tool.diagnostic as diagnostic
import geneset_network_analysis_tool.painter as paint
from geneset_network_analysis_tool.utils import Converter
import pandas as pd
import scipy
import time
import scipy.linalg.interpolative
from copy import copy, deepcopy
import itertools
import tables
import seaborn as sns
from matplotlib import pyplot as plt


def fancy_dendrogram(*args, **kwargs):
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        plt.title('Hierarchical Clustering Dendrogram (truncated)')
        plt.xlabel('sample index or (cluster size)')
        plt.ylabel('distance')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')
    return ddata

def __load_network(filename):
    """ Loads network from file
    """
    if filename.endswith('tsv'):
        network = ps.TSVParser().read(filename)
    elif filename.endswith('gpickle'):
        network = nx.read_gpickle(filename)
    else:
        network = ps.Tab2Parser().read(filename)

    return network


def __load_geneset(filename, setname=None):
    """Loads a geneset from file
    """

    geneset = ps.GMTParser().read(filename)
    if setname:
        if setname in geneset:
            temp = geneset[setname]
            geneset.clear()
            geneset[setname] = temp
        else:
            logging.error("Cannot find geneset: %s" % setname)
            sys.exit(-1)

    return geneset


# def __convert_entrez2symbol(entrez_data):
#    with open("/primary_data/Homo_sapiens.gene_info","r") as f:
#        table_id=pd.read_table(f, usecols=["GeneID", "Symbol"])
#    table_id=table_id.set_index("GeneID")
#    symbol_data=entrez_data.map(table_id.Symbol)
#    return symbol_data

def __read_RW_matrix(RW_matrix_filename):
    """ Reads matrix from hdf5 file
    """

    if RW_matrix_filename.endswith('hdf5'):
        # hdf5_file="/home/viola/Documents/UniEssex/repos/geneset-network-analysis/processed_data/barabasi_distance_matrix.hdf5"
        hdf5_nodes, hdf5_data = ps.Hdf5MatrixParser().read(RW_matrix_filename)
    else:
        logging.error("invalid input format for matrix")

    return hdf5_nodes, hdf5_data


def __read_distance_matrix(distance_matrix_filename):
    """ Reads the large matrix .lm text file for the diffusion matrix and creates
    a hdf5 file to work with it
    """
    # filename="/home/viola/Documents/UniEssex/repos/geneset-network-analysis/processed_data/barabasi_distance_matrix.lm"
    if distance_matrix_filename.endswith('lm'):

        # with open(distance_matrix_filename) as f:
        #    nodes_names = f.readline()[:-7].replace("'","").replace("\t",",").replace(" ","").split(",")
        # hdf5_file.close()
        hdf5_path = "/home/viola/Documents/UniEssex/repos/geneset-network-analysis/processed_data/barabasi_distance_matrix.hdf5"
        with tables.open_file(hdf5_path, mode='w') as hdf5_file:

            line_index = 0
            with open(distance_matrix_filename, "r") as f:
                for line in f:
                    if line_index == 0:
                        nodes_name = line[:-7].replace("'", "").replace(
                            "\t", ",").replace(" ", "").split(",")
                        # TODO ricalcola matrix e rimuovi -7
                        hdf5_nodes = hdf5_file.create_array(
                            hdf5_file.root, 'nodes', nodes_name)

                        line_index += 1
                        hdf5_data = hdf5_file.create_array(
                            hdf5_file.root, 'matrix', np.zeros((len(nodes_name), len(nodes_name))))
                    else:
                        i, j, sp = [int(x)
                                    for x in line.replace("\n", "").split("\t")]
                        hdf5_data[i, j] = sp
                        line_index += 1
                        if line_index % 13329 == 0:
                            print(line_index // 13329)

    elif distance_matrix_filename.endswith('hdf5'):
        # hdf5_file="/home/viola/Documents/UniEssex/repos/geneset-network-analysis/processed_data/barabasi_distance_matrix.hdf5"
        hdf5_nodes, hdf5_data = ps.Hdf5MatrixParser().read(distance_matrix_filename)
    else:
        logging.error("invalid input format for matrix")

    return hdf5_nodes, hdf5_data


def network_summary(network_file, output_folder, output_name, geneset_file=None, setname=None):
    """
    This function saves the principal info on a graph:
    - network properties
    - degree distribution
    """

    network = __load_network(network_file)
    D=dict(nx.degree(network))

    # TODO: organise code better
    if geneset_file:
        geneset=__load_geneset(geneset_file,setname)
        for setname, item in geneset.items():


            graph = nx.subgraph(network, item)

            degree = np.array(list(dict(nx.degree(graph)).values()))

            n_nodes = nx.number_of_nodes(graph)
            n_edges = nx.number_of_edges(graph)

            if n_nodes > 1:

                density = (2 * n_edges) / ((n_nodes) * (n_nodes - 1))

            else: density=0

            degrees={k: v for k, v in D.items() if k in item}
            degrees = sorted(degrees.items(), key=lambda kv: kv[1])




            with open(output_folder +output_name + setname + "graph_summary.txt", "w") as file1:
                file1.write("Info: " + nx.info(graph))
                file1.write("\nDensity: " + str(density))
                file1.write("\nTop degrees: " + str(degrees[-10:]))
                file1.write("\nmin degree = " + str(np.min(degree)))
                file1.write("\nmax degree = " + str(np.max(degree)))
                file1.write("\nmedian degree = " + str(np.median(degree)))
                file1.write("\ndegree mode = " + str(scipy.stats.mode(degree)))
                file1.write("\ndisconnected nodes = " + str(np.sum(degree == 0)))
                file1.write("\nClustering: \naverage clsutering" +
                            str(nx.average_clustering(graph)))

    else:

        graph = network

        degree = np.array(list(dict(nx.degree(graph)).values()))

        n_nodes = nx.number_of_nodes(graph)
        n_edges = nx.number_of_edges(graph)

        degrees={k: v for k, v in D.items()}
        degrees = sorted(degrees.items(), key=lambda kv: kv[1])

        density = (2 * n_edges) / ((n_nodes) * (n_nodes - 1))

        with open(output_folder +output_name + "_graph_summary.txt", "w") as file1:
            file1.write("Info: " + nx.info(graph))
            file1.write("\nDensity: " + str(density))

            file1.write("\nmin degree = " + str(np.min(degree)))
            file1.write("\nmax degree = " + str(np.max(degree)))
            file1.write("\nmedian degree = " + str(np.median(degree)))
            file1.write("\ndegree mode = " + str(scipy.stats.mode(degree)))
            file1.write("\ndisconnected nodes = " + str(np.sum(degree == 0)))
            file1.write("\nClustering: \naverage clsutering" +
                            str(nx.average_clustering(graph)))

        f, axis = plt.subplots(1)
        sns.distplot(degree, ax=axis)
        f.savefig(output_folder + output_name
         + "_degree_distribution.pdf", format='pdf')


def summary_composition():
    pass

################################################################################
######### SINGLE SET ANALYSES ##################################################
################################################################################
def analyse_total_degree( network_file: "network file",
                            geneset_file: "GMT geneset file",
                            outpath: "output folder where to place all generated files",
                            prefix: "prefix to add to all generated files",
                            setname: "Geneset to analyse" = None,
                            number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
                            cores: "Number of cores for the multiprocessing" = 1,
                            show_results: "barplot of results" = False,
                            show_null: "plot null distribution" = False,
                            symbol: "True if we want to print th output in symbol names" = False):

    '''
        Performs the analysis of internal degree.

        It computes a p-value for the ratio of internal degree
        of the geneset being bigger than the one expected by chance
        for a geneset of the same size.
    '''

    network = __load_network(network_file)

    geneset = __load_geneset(geneset_file, setname)

    setnames = [key for key in geneset.keys()]

    output1 = out.Output(network_file, outpath, prefix,
                         "analyse_total_degree", geneset_file, setnames)
    logging.info("Output-folder= " + output1.output)
    output1.create_st_table_empirical("table_total_degree")

    st_test = st.StatisticalTest(st.geneset_total_degree_statistic, network)

    for setname, item in geneset.items():

        start = time.time()
        item = set(item)

        observed, pvalue, null_d, n_mapped, n_geneset = st_test.empirical_pvalue(
            item, max_iter=number_of_permutations, alternative='greater', cores= cores)
        # observed_z=(observed-np.mean(null_d))/np.std(null_d)
        end = time.time()

        logging.info("Setname:" + setname)
        logging.info("Observed: %g p-value: %g" % (observed, pvalue))
        logging.info("Null mean: %g null variance: %g" %
                     (np.mean(null_d), np.var(null_d)))
        logging.info("Time: %g" % (end - start))

        output1.update_st_table_empirical(
            setname, n_mapped, n_geneset, number_of_permutations, observed, pvalue, np.mean(null_d), np.var(null_d))
        if show_null:
            diagnostic.plot_null_distribution(
                null_d, observed, output1.output, setname=setname)

    output1.save_output_summary()

    if show_results:
        paint.paint_datasets_stats(
            output1.output_table, output1.output, "total_degree")

def analyse_internal_degree( network_file: "network file",
                            geneset_file: "GMT geneset file",
                            outpath: "output folder where to place all generated files",
                            prefix: "prefix to add to all generated files",
                            setname: "Geneset to analyse" = None,
                            number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
                            cores: "Number of cores for the multiprocessing" = 1,
                            show_results: "barplot of results" = False,
                            show_null: "plot null distribution" = False,
                            symbol: "True if we want to print th output in symbol names" = False):

    '''
        Performs the analysis of internal degree.

        It computes a p-value for the ratio of internal degree
        of the geneset being bigger than the one expected by chance
        for a geneset of the same size.
    '''

    network = __load_network(network_file)

    geneset = __load_geneset(geneset_file, setname)

    setnames = [key for key in geneset.keys()]

    output1 = out.Output(network_file, outpath, prefix,
                         "analyse_internal_degree", geneset_file, setnames)
    logging.info("Output-folder= " + output1.output)
    output1.create_st_table_empirical("table_internal_degree")

    st_test = st.StatisticalTest(st.geneset_internal_degree_statistic, network)

    for setname, item in geneset.items():

        start = time.time()
        item = set(item)

        observed, pvalue, null_d, n_mapped, n_geneset = st_test.empirical_pvalue(
            item, max_iter=number_of_permutations, alternative='greater', cores= cores)
        # observed_z=(observed-np.mean(null_d))/np.std(null_d)
        end = time.time()

        logging.info("Setname:" + setname)
        logging.info("Observed: %g p-value: %g" % (observed, pvalue))
        logging.info("Null mean: %g null variance: %g" %
                     (np.mean(null_d), np.var(null_d)))
        logging.info("Time: %g" % (end - start))

        output1.update_st_table_empirical(
            setname, n_mapped, n_geneset, number_of_permutations, observed, pvalue, np.mean(null_d), np.var(null_d))
        if show_null:
            diagnostic.plot_null_distribution(
                null_d, observed, output1.output, setname=setname)

    output1.save_output_summary()

    if show_results:
        paint.paint_datasets_stats(
            output1.output_table, output1.output, "internal_degree")

def analyse_RW(network_file: "network file, use a network with weights",
               geneset_file: "GMT geneset file",
               RW_dict_file: "pickle obtained calculating a random walk matrix",
               outpath: "output folder where to place all generated files",
               prefix: "prefix to add to all generated files",
               setname: "Geneset to analyse" = None,
               weight: "RW" = "RW",
               number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
               cores: "Number of cores for the multiprocessing" = 1,
               show_matrix: "plotting flag, if true the diffusion matrix for each geneset is saved " = False,
               show_results: "barplot of results" = False,
               show_null: "plot null distribution" = False):
    """
        Performs the analysis of random walk probabilities. Given the RW matrix ( either normal random walk or RW with restart),
        it compares the probability of walking between thhe genes in the geneset compared to those of walking between the nodes
        of a geneset with the same size
    """

    network = __load_network(network_file)
    network = nx.Graph(network.subgraph(
        max(nx.connected_components(network), key=len)))

    geneset = __load_geneset(geneset_file, setname)
    RW_dict = {} #ps.DiffusionDictParser().read(RW_dict_file)

    RW_dict["nodes"], RW_dict["matrix"] = __read_distance_matrix(RW_dict_file)

    setnames = [key for key in geneset.keys()]
    output1 = out.Output(network_file, outpath, prefix,
                         "analyse_RW", geneset_file, setnames)

    logging.info("Output-folder= " + output1.output)
    output1.create_st_table_empirical("table_RW")

    output1.set_diffusion_matrix(RW_dict_file)

    st_test = st.StatisticalTest(st.geneset_RW_statistic, network, RW_dict)

    for setname, item in geneset.items():

        item = set(item).intersection(set(list(network.nodes)))
        observed, pvalue, null_d, n_mapped, n_geneset = st_test.empirical_pvalue(
            item, max_iter=number_of_permutations, alternative='greater')
        # observed_z=(observed-np.mean(null_d))/np.std(null_d
        logging.info("Plotting diagnostic"
                     + str(output1.output))

        if (len(item)>0 and show_null):
            diagnostic.plot_null_distribution(
                null_d, observed, output1.output, setname=setname)

        logging.info("Setname:" + setname)
        logging.info("Observed: %g p-value: %g" % (observed, pvalue))
        logging.info("Null mean: %g null variance: %g" %
                     (np.mean(null_d), np.var(null_d)))
        output1.update_st_table_empirical(
            setname, n_mapped, n_geneset, number_of_permutations, observed, pvalue, np.mean(null_d), np.var(null_d))

        logging.info("show matrix: "+str(show_matrix) )
        if (len(item)>0 and show_matrix):
            print(len(item))
            print(item)
            painter_rw = paint.Painter_RW(
                network, output1.output , setname, RW_dict, item)
            logging.info("Painting matrix")
            painter_rw.plot_matrix()

    output1.save_output_summary()

    if show_results:
        paint.paint_datasets_stats(
            output1.output_table, output1.output, "RW")

def analyse_module(network_file: "network file",
                   geneset_file: "GMT geneset file",
                   outpath: "output folder where to place all generated files",
                   prefix: "prefix to add to all generated files",
                   setname: "Geneset to analyse" = None,
                   number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
                   cores: "Number of cores for the multiprocessing" = 1,
                   create_output_LCC: "flag for creating a GMT file with the LCC lists" = False,
                   show_results: "barplot of results" = False,
                   show_null: "plot null distribution" = False,
                   symbol: "True if we want to print th output in symbol names" = False):
    '''
        Performs network module analysis.

        It computes a p-value for the largest connected component
        of the geneset being bigger than the one expected by chance
        for a geneset of the same size.
    '''
    network = __load_network(network_file)

    geneset = __load_geneset(geneset_file, setname)

    setnames = [key for key in geneset.keys()]

    output1 = out.Output(network_file, outpath, prefix,
                         "analyse_module", geneset_file, setnames)
    logging.info("Output-folder= " + output1.output)
    output1.create_st_table_empirical("table_module")

    st_test = st.StatisticalTest(st.geneset_module_statistic, network)

    for setname, item in geneset.items():

        start = time.time()
        item = set(item)
        if create_output_LCC:

            module = nx.subgraph(network, item)
            if len(module.nodes)>0:
                LCC = sorted(list(nx.connected_components(module)),
                            key=len, reverse=True)[0]
            else:
                LCC=[]
            if symbol:
                converter = Converter()
                LCC = converter.entrez2symbol(LCC)

            output1.add_GMT_entry(setname, "module_analysis", LCC)

        observed, pvalue, null_d, n_mapped, n_geneset = st_test.empirical_pvalue(
            item, max_iter=number_of_permutations, alternative='greater')
        # observed_z=(observed-np.mean(null_d))/np.std(null_d)
        end = time.time()

        logging.info("Setname:" + setname)
        logging.info("Observed: %g p-value: %g" % (observed, pvalue))
        logging.info("Null mean: %g null variance: %g" %
                     (np.mean(null_d), np.var(null_d)))
        logging.info("Time: %g" % (end - start))

        output1.update_st_table_empirical(
            setname, n_mapped, n_geneset, number_of_permutations, observed, pvalue, np.mean(null_d), np.var(null_d))
        if show_null:
            diagnostic.plot_null_distribution(
                null_d, observed, output1.output, setname=setname)

    output1.save_output_summary()

    if create_output_LCC:
        output1.create_GMT_output()

    if show_results:
        paint.paint_datasets_stats(
            output1.output_table, output1.output, "module")

def analyse_location(network_file: "network file",
                     distance_matrix_filename: "distance matrix file generated by geneset_network_analysis_tool",
                     geneset_file: "GMT geneset file",
                     outpath: "output folder where to place all generated files",
                     prefix: "prefix to add to all generated files",
                     setname: "Geneset to analyse" = None,
                     number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
                     cores: "Number of cores for the multiprocessing" = 1,
                     show_matrix: "plotting flag, if true the distance matrix for each geneset is saved " = False,
                     show_results: "barplot of results" = False):
    '''
        Performs network location analysis.

        It computes a p-value for the average shortest path length
        of the geneset being smaller than expected by chance
        for a geneset of the same size.
    '''

    network = __load_network(network_file)
    network = nx.Graph(network.subgraph(
        max(nx.connected_components(network), key=len)))
    geneset = __load_geneset(geneset_file, setname)
    #distance_matrix = ps.DistanceMatrixParser().read(distance_matrix_file)
    diz = {}
    diz["nodes"], diz["matrix"] = __read_distance_matrix(
        distance_matrix_filename)

    setnames = [key for key in geneset.keys()]

    output1 = out.Output(distance_matrix_filename, outpath, prefix,
                         "analyse_location", geneset_file, setnames)
    logging.info("Output-folder= " + output1.output)
    output1.create_st_table_empirical("table_SP")

    st_test = st.StatisticalTest(
        st.geneset_localisation_statistic, network, diz)

    for setname, item in geneset.items():

        item = set(item).intersection(set(list(network.nodes)))
        observed, pvalue, null_d, n_mapped, n_geneset = st_test.empirical_pvalue(
            item, max_iter=number_of_permutations)

        logging.info("Observed: %g p-value: %g" % (observed, pvalue))
        logging.info("Null mean: %g null variance: %g" %
                     (np.mean(null_d), np.var(null_d)))
        logging.info("1th percentile: %g " % (np.percentile(null_d, 1)))

        output1.update_st_table_empirical(
            setname, n_mapped, n_geneset, number_of_permutations, observed, pvalue, np.mean(null_d), np.var(null_d))

    output1.save_output_summary()

    if show_results:
        paint.paint_datasets_stats(
            output1.output_table, output1.output, "SP")

def analyse_location_median(network_file: "network file",
                            distance_matrix_filename: "distance matrix file generated by geneset_network_analysis_tool",
                            geneset_file: "GMT geneset file",
                            outpath: "output folder where to place all generated files",
                            prefix: "prefix to add to all generated files",
                            setname: "Geneset to analyse" = None,
                            number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
                            cores: "Number of cores for the multiprocessing" = 1,
                            show_matrix: "plotting flag, if true the distance matrix for each geneset is saved " = False,
                            show_results: "barplot of results" = False):
    '''
        Performs network location analysis.

        It computes a p-value for the average shortest path length
        of the geneset being smaller than expected by chance
        for a geneset of the same size.
    '''

    network = __load_network(network_file)
    network = nx.Graph(network.subgraph(
        max(nx.connected_components(network), key=len)))
    geneset = __load_geneset(geneset_file, setname)
    #distance_matrix = ps.DistanceMatrixParser().read(distance_matrix_file)
    diz = {}
    diz["nodes"], diz["sp_matrix"] = __read_distance_matrix(
        distance_matrix_filename)

    setnames = [key for key in geneset.keys()]

    output1 = out.Output(distance_matrix_filename, outpath, prefix,
                         "analyse_location_median", geneset_file, setnames)
    logging.info("Output-folder= " + output1.output)
    output1.create_st_table_empirical("table_SP_median")

    st_test = st.StatisticalTest(
        st.geneset_localisation_statistic_median, network, diz)

    for setname, item in geneset.items():

        item = set(item).intersection(set(list(network.nodes)))
        observed, pvalue, null_d, n_mapped, n_geneset = st_test.empirical_pvalue(
            item, max_iter=number_of_permutations)

        logging.info("Observed: %g p-value: %g" % (observed, pvalue))
        logging.info("Null mean: %g null variance: %g" %
                     (np.mean(null_d), np.var(null_d)))
        logging.info("1th percentile: %g " % (np.percentile(null_d, 1)))

        output1.update_st_table_empirical(
            setname, n_mapped, n_geneset, number_of_permutations, observed, pvalue, np.mean(null_d), np.var(null_d))

    output1.save_output_summary()

    if show_results:
        paint.paint_datasets_stats(
            output1.output_table, output1.output, "SP_median")

################################################################################
######### Degree Distribution test #############################################
################################################################################

def test_degree_distribution(network_file: "network file",
                   geneset_file: "GMT geneset file",
                   outpath: "output folder where to place all generated files",
                   prefix: "prefix to add to all generated files",
                   setname: "Geneset to analyse" = None,
                   show_results: "barplot of results" = False):
    '''
        Performs degree distribution test
    '''
    network = __load_network(network_file)

    geneset = __load_geneset(geneset_file, setname)

    setnames = [key for key in geneset.keys()]

    output1 = out.Output(network_file, outpath, prefix,
                         "test_degree", geneset_file, setnames)
    logging.info("Output-folder= " + output1.output)
    output1.create_st_table_empirical("table_degree")

    st_test = KS.KSTest(KS.degree_distribution, network)

    for setname, item in geneset.items():

        item = set(item)

        observed, pvalue, n_mapped, n_geneset = st_test.apply_test(item)
        # observed_z=(observed-np.mean(null_d))/np.std(null_d)

        logging.info("Setname:" + setname)
        logging.info("Observed: %g p-value: %g" % (observed, pvalue))

        output1.update_st_table_empirical(
            setname, n_mapped, n_geneset, 1, observed, pvalue, 0, 0)

    output1.save_output_summary()


    if show_results:
        paint.paint_datasets_stats(
            output1.output_table, output1.output, "test_degree")

################################################################################
######### COMPARISONS ##########################################################
################################################################################


def comparison_shortest_path(network_file: "network file",
                             distance_matrix_filename: "distance matrix file generated by geneset_network_analysis_tool",
                             A_geneset_file: "GMT geneset file, if it's the only parameter passed the analysis is gonna be run on all the couples of datasets, otherwise specify the other files and setnames",
                             outpath: "output folder where to place all generated files",
                             prefix: "prefix to add to all generated files",
                             setname_A: "Geneset A to analyse" = None,
                             B_geneset_file: "GMT geneset file" = None,
                             setname_B: "Geneset B to analyse" = None,
                             cores: "Number of cores for the multiprocessing" = 1,
                             number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
                             show_matrix: "plotting flag, if true the diffusion matrix for each pair of genesets is saved" = False,
                             show_results: "heatmap of results" = False):
    '''
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
    '''

    network = __load_network(network_file)
    network = nx.Graph(network.subgraph(
        max(nx.connected_components(network), key=len)))

    sp_diz = {}
    sp_diz["nodes"], sp_diz["matrix"] = __read_distance_matrix(
        distance_matrix_filename)

    if (setname_A and setname_B == None and B_geneset_file == None):
        logging.error(" this analysis requires at least two genesets ")

    geneset_A = __load_geneset(A_geneset_file, setname_A)

    if B_geneset_file:
        geneset_B = __load_geneset(B_geneset_file, setname_B)
    else:
        if setname_B:
            geneset_B = __load_geneset(A_geneset_file, setname_B)
        else:
            geneset_B = None

    st_comparison = sc.StatisticalComparison(
        sc.comparison_shortest_path, network, diz=sp_diz, n_proc=cores)

    if not geneset_B:  # Analysis of genesets inside a single file

        logging.info("Analysising all the sets in " + A_geneset_file)
        setnames = [key for key in geneset_A.keys()]
        logging.info("Setnames: " + str(setnames))

        # Creating the output table
        output1 = out.Output(
            network_file, outpath, prefix, "comparison_shortest_path", A_geneset_file, setnames)
        output1.add_output_text(
            " distance matrix = " + str(distance_matrix_filename))
        logging.info("Output-folder= " + output1.output)
        output1.create_comparison_table_empirical(
            "table_comparison_shortest_path")

        for pair in itertools.combinations(setnames, 2):
            logging.info("Analysing " + str(pair[0]) + " and " + str(pair[1]))
            overlaps=set(geneset_A[pair[0]]).intersection(set(geneset_A[pair[1]]))
            logging.info('There are %d genes shared between A and B' %len(overlaps))
            observed, pvalue, null_d, A_mapped, B_mapped = st_comparison.comparison_empirical_pvalue(
                set(geneset_A[pair[0]]), set(geneset_A[pair[1]]), max_iter=number_of_permutations)
            # Save the results
            output1.update_comparison_table_empirical(pair[0], pair[1],
                                                      len(set(
                                                          geneset_A[pair[0]])), A_mapped,
                                                      len(set(
                                                          geneset_A[pair[1]])), B_mapped, len(overlaps),
                                                      number_of_permutations, observed, pvalue,
                                                      np.mean(null_d), np.var(null_d))

    else:  # Analysis of genesets into two different gmt files

        logging.info("geneset_A contains %d sets", (len(geneset_A)))
        sets_A = [key for key in geneset_A.keys()]
        logging.info("Setnames in A: " + str(sets_A))
        logging.info("geneset_B contains %d sets", (len(geneset_B)))
        sets_B = [key for key in geneset_B.keys()]
        logging.info("Setnames in B: " + str(sets_B))

        output1 = out.Output(network_file, outpath, prefix,  "comparison_shortest_path",
                             A_geneset_file, sets_A, B_geneset_file, sets_B)
        output1.add_output_text(
            " distance matrix = " + str(distance_matrix_filename))
        logging.info("Output-folder= " + output1.output)
        output1.create_comparison_table_empirical(
            "table_comparison_shortest_path")

        for set_A, item_A in geneset_A.items():
            for set_B, item_B in geneset_B.items():
                n_overlaps = len(set(item_A).intersection(set(item_B)))

                observed, pvalue, null_d, A_mapped, B_mapped = st_comparison.comparison_empirical_pvalue(
                    set(item_A), set(item_B), max_iter=number_of_permutations)

                logging.info("Observed: %g p-value: %g" % (observed, pvalue))
                logging.info("Null mean: %g null variance: %g" %
                             (np.mean(null_d), np.var(null_d)))
                logging.info("1th percentile: %g " %
                             (np.percentile(null_d, 1)))

                output1.update_comparison_table_empirical(set_A, set_B,
                                                          len(set(item_A)
                                                              ), A_mapped,
                                                          len(set(item_B)
                                                              ), B_mapped,
                                                          n_overlaps,
                                                          number_of_permutations, observed, pvalue,
                                                          np.mean(null_d), np.var(null_d))

    output1.save_output_summary()  # Save Summary of Analysis

    if show_results:
        paint.paint_comparison_stats(
            output1.output_table, output1.output + "/", "sp")


def comparison_random_walk(network_file: "network file",
                           A_geneset_file: "GMT geneset file",
                           RW_dict_file: "pickle obtained calculating a random walk matrix",
                           outpath: "output folder where to place all generated files",
                           prefix: "prefix to add to all generated files",
                           setname_A: "Geneset A to analyse" = None,
                           B_geneset_file: "GMT geneset file" = None,
                           setname_B: "Geneset B to analyse" = None,
                           cores: "Number of cores for the multiprocessing" = 1,
                           number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
                           show_matrix: "plotting flag, if true the diffusion matrix for each pair of genesets is saved " = False,
                           show_results: "heatmap of results" = False):
    '''
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
    '''

    network = __load_network(network_file)
    network = nx.Graph(network.subgraph(
        max(nx.connected_components(network), key=len)))

    # Read datasets
    if (setname_A and setname_B == None and B_geneset_file == None):
        logging.error(" this analysis requires at least two genesets ")

    logging.info("Reading diffusion matrix")
    RW_dict = {}
    RW_dict["nodes"], RW_dict["matrix"] = __read_distance_matrix(
        RW_dict_file)

    geneset_A = __load_geneset(A_geneset_file, setname_A)

    if B_geneset_file:
        geneset_B = __load_geneset(B_geneset_file, setname_B)
    else:
        if setname_B:
            geneset_B = __load_geneset(A_geneset_file, setname_B)
        else:
            geneset_B = None

    # define analysis and network  (network, genesetA, genesetB, diz={})
    st_comparison = sc.StatisticalComparison(
        sc.comparison_random_walk, network, n_proc=cores, diz=RW_dict)

    if not geneset_B:

        logging.info("Analysising all the sets in " + A_geneset_file)
        setnames = [key for key in geneset_A.keys()]
        logging.info("Setnames: " + str(setnames))

        output1 = out.Output(network_file, outpath, prefix,
                             "comparison_random_walk", A_geneset_file, setnames)
        output1.add_output_text(" RW matrix = " + str(RW_dict_file))
        logging.info("Output-folder= " + output1.output)
        output1.create_comparison_table_empirical(
            "table_comparison_random_walk")

        for pair in itertools.combinations(setnames, 2):
            logging.info("Analysing " + str(pair[0]) + " and " + str(pair[1]))

            n_overlaps = len(set(geneset_A[pair[0]]).intersection(
                set(geneset_A[pair[1]])))

            observed, pvalue, null_d, A_mapped, B_mapped = st_comparison.comparison_empirical_pvalue(set(
                geneset_A[pair[0]]), set(geneset_A[pair[1]]), max_iter=number_of_permutations, alternative='greater')

            output1.update_comparison_table_empirical(pair[0], pair[1],
                                                      len(set(
                                                          geneset_A[pair[0]])), A_mapped,
                                                      len(set(
                                                          geneset_A[pair[1]])), B_mapped,
                                                      n_overlaps,
                                                      number_of_permutations, observed, pvalue,
                                                      np.mean(null_d), np.var(null_d))

    else:

        logging.info("geneset_A contains %d sets" % ((len(geneset_A))))
        sets_A = [key for key in geneset_A.keys()]
        logging.info("Setnames in A: " + str(sets_A))
        logging.info("geneset_B contains %d sets" % ((len(geneset_B))))
        sets_B = [key for key in geneset_B.keys()]
        logging.info("Setnames in B: " + str(sets_B))

        output1 = out.Output(network_file, outpath, prefix,"comparison_random_walk",
                             A_geneset_file, sets_A, B_geneset_file, sets_B)
        logging.info("Output-folder= " + output1.output)
        output1.create_comparison_table_empirical("table_random_walk")
        output1.add_output_text(" RW matrix = " + str(RW_dict_file))

        for set_A, item_A in geneset_A.items():
            for set_B, item_B in geneset_B.items():
                n_overlaps = len(set(item_A).intersection(set(item_B)))

                observed, pvalue, null_d, A_mapped, B_mapped = st_comparison.comparison_empirical_pvalue(
                    set(item_A), set(item_B), max_iter=number_of_permutations)

                logging.info("Observed: %g p-value: %g" % (observed, pvalue))
                logging.info("Null mean: %g null variance: %g" %
                             (np.mean(null_d), np.var(null_d)))
                logging.info("1th percentile: %g " %
                             (np.percentile(null_d, 1)))

                output1.update_comparison_table_empirical(set_A, set_B,
                                                          len(set(item_A)
                                                              ), A_mapped,
                                                          len(set(item_B)
                                                              ), B_mapped,
                                                          n_overlaps,
                                                          number_of_permutations, observed, pvalue,
                                                          np.mean(null_d), np.var(null_d))

    output1.save_output_summary()  # Save Summary of Analysis
    if show_results:
        paint.paint_comparison_stats(
            output1.output_table, output1.output, "RWR")

################################################################################
######### BUILDING FUNCTIONS ###################################################
################################################################################


def build_distance_matrix(network_file: "network file",
                          output_file: "distance matrix output file, use a format between .lm.txt and .hdf5",
                          giant_component_only: "compute the shortest paths only for nodes in the giant component" = True
                          ):
    '''
        Build a shortest path distance matrix for a given network.
        Matrix can be saved as a lm.txt file or a .hdf5 one.
    '''

    #network = __load_network(network_file)
    network = __load_network(network_file)
    if giant_component_only:
        network = network.subgraph(max(nx.connected_components(network), key=len))

    print(type(network))
    distance_matrix = nx.all_pairs_shortest_path_length(network)
    print(type(distance_matrix))

    nodes = list(network.nodes())

    if output_file.endswith(".hdf5"):
        with tables.open_file(output_file, mode='w') as hdf5_file:

            # create a hdf5 file with two objects:
            # - one is the nodes array,
            hdf5_nodes = hdf5_file.create_array(hdf5_file.root, 'nodes', nodes)
            # -  the other is the shorthest path distance matrix
            hdf5_data = hdf5_file.create_array(
                hdf5_file.root, 'matrix', np.zeros((len(nodes), len(nodes))))

            for node_i, k in distance_matrix:
                i = nodes.index(node_i)

                for node_j, sp in k.items():
                    j = nodes.index(node_j)
                    if j >= i: #saves only the upper triangular matrix
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

    #pickle.dump(distance_matrix, open(output_file, 'wb'))


def build_RWR_diffusion(network_file: "network file",
                        beta = 0.85,
                        output_file: "distance matrix output file (use .hdf5) " = None,
                        output_figure: "diffusion matrix output filename" = None):
    '''
        Build the RWR_diffusion_matrix
    '''

    network = __load_network(network_file)
    network = network.subgraph(max(nx.connected_components(network), key=len))

    nodes = list(network.nodes())

    logging.info("Beginning to calculate RWR matrix")

    A = nx.adjacency_matrix(network)
    K = 1 / np.array(list(dict(network.degree()).values()))
    D = scipy.sparse.dia_matrix((K, [0]), shape=(A.shape))
    A = A.dot(D)
    n = np.shape(A)[1]

    if output_file.endswith(".hdf5"):
        with tables.open_file(output_file, mode='w') as hdf5_file:

            # create a hdf5 file with two objects:
            # - one is the nodes array,
            hdf5_nodes = hdf5_file.create_array(hdf5_file.root, 'nodes', nodes)
            # -  the other is the shorthest path distance matrix
            hdf5_data = hdf5_file.create_array(
                hdf5_file.root, 'matrix', beta * np.linalg.inv(np.eye(n) - (1. - beta) * A) )

            logging.info("Saving network")
    else:
        return (beta * np.linalg.inv(np.eye(n) - (1. - beta) * A) )

    if output_figure:
        diagnostic.plot_diffusion_matrix(nodes, hdf5_data, output_figure)

def build_random_walk(network_file: "network file",
                      output_file: "distance matrix output file (pickle format)" = None,
                      output_figure: "diffusion matrix output filename" = None):
    '''
        Build the random walk diffusion_matrix
    '''

    network = __load_network(network_file)
    network = network.subgraph(max(nx.connected_components(network), key=len))

    output_dict = {}
    output_dict["nodes"] = list(network.nodes())

    logging.info("Creating the RW network")

    D = np.array([float(i) for i in dict(network.degree()).values()])
    D = D / np.sum(D)
    print(D)
    output_dict["RW"] = np.repeat(D.T[:, np.newaxis], len(D), axis=1)

    logging.info("Saving network")

    if output_figure:
        diagnostic.plot_diffusion_matrix(
            output_dict["nodes"], output_dict["RW"], output_figure)

    if output_file:
        pickle.dump(output_dict, open(output_file + ".pickle", 'wb'))
    else:
        return output_dict


def build_hitting_time_matrix(network_file: "network file",
                              output_file: "distance matrix output file (pickle format)" = None):
    '''
        Build the hitting time matrix
    '''

    network = __load_network(network_file)
    network = network.subgraph(max(nx.connected_components(network), key=len))

    output_dict = {}
    output_dict["nodes"] = network.nodes()
    output_dict["tau"] = 0

    A = nx.adjacency_matrix(network)
    print(type(A))
    K = 1 / np.array(list(dict(network.degree()).values()))
    D = scipy.sparse.dia_matrix((K, [0]), shape=(A.shape))
    A = A.dot(D)
    print(type(A))
    print("Shape of A: " + str(A.shape))

    logging.info("Beginning to calculate hitting time matrix")
    output_dict["tau"] = np.zeros((A.shape))

    for i in np.arange(A.shape[0]):

        print("node=" + str(i))

        A1 = A.copy().tocsc()
        A1[:, i] = 0
        A1[i, i] = 1
        A1 = A1.tocsr()
        Wk = A1[i, :].copy()

        H = np.zeros((1, A1.shape[0]))
        F = Wk

        diff = 10

        k = 0
        while np.max(diff) > 0.0001:
            # print(k)
            temp = deepcopy(Wk)
            Wk = np.dot(Wk, A1)
            diff = Wk - temp
            F = F + (k + 2) * diff
            k += 1

        output_dict["tau"][:, i] = F.toarray()

    logging.info("Hitting time matrix calculated")
    if output_file:
        pickle.dump(output_dict, open(output_file, 'wb'))
    else:
        return tau


def build_graph(network_file: "network file",
                geneset_file: "geneset file",
                output_folder: "graphml network for visualisation",
                setname: "setname" = None,
                giant_component_only: "compute the shortest paths only for nodes in the giant component" = True):
    '''
        Build a shortest path distance matrix for a given network.
    '''
    network = __load_network(network_file)
    geneset = __load_geneset(geneset_file, setname)

    if giant_component_only:
        network = network.subgraph(
            max(nx.connected_components(network), key=len))

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

        nx.write_graphml(new_network_minimal, output_folder
                         + setname + "_minimal.graphml")


################################################################################
######### CLUSTERING ##########################################################
################################################################################


def hierarchical_clustering_RW(network_file: "network file, use a network with weights",
                                geneset_file: "GMT geneset file",
                                RW_dict_file: "pickle obtained calculating a random walk matrix",
                                outpath: "output folder where to place all generated files",
                                prefix: "prefix to add to all generated files",
                                setname: "Geneset to analyse" = None,
                                weight: "RW" = "RW",
                                show_results: "image of the dendrogram" = False,
                                show_labels=False):

    """
        The script performs the hierarchical clustering of the RW matrix of a geneset
    """

    network = __load_network(network_file)
    network = nx.Graph(network.subgraph(
          max(nx.connected_components(network), key=len)))

    geneset = __load_geneset(geneset_file, setname)

    RW_dict = {}
    RW_dict["nodes"], RW_dict["matrix"] = __read_distance_matrix(RW_dict_file)

    #RW_dict = ps.DiffusionDictParser().read(RW_dict_file)

    logging.info(RW_dict.keys())

    setnames = [key for key in geneset.keys()]

    #create output
    output1 = out.Output(network_file, outpath, prefix,
                         "hierarchical clustering", geneset_file, setnames)
    logging.info("Output-folder= " + output1.output)
    output1.set_diffusion_matrix(RW_dict_file)

    # Clustering for each dataset
    for setname, item in geneset.items():
        logging.info("\nSetname: " + setname)
        print(item)
        print(RW_dict["nodes"])

        #item = set(item.encode('itf-8')).intersection(set(list(RW_dict["nodes"])))

        #mapping each gene to the ones of the matrix
        #geneset_index=[RW_dict["nodes"].index(i.encode('utf-8')) for i in item if i.encode('utf-8') in RW_dict["nodes"]]
        geneset_index=[RW_dict["nodes"].index(i)  for i in RW_dict["nodes"]]
        item = [i for i in RW_dict["nodes"]]

        # Allocate a matrix, it is going to be the "distance" matrix expressed as 1-Probability
        X = np.ones((len(geneset_index),len(geneset_index)))
        print(X)
        for i in range(len(geneset_index)):
            for j in range(len(geneset_index)):
                if i!=j:
                    # we use the maximum probability of interaction between i and j
                    X[i,j]=1-max(RW_dict["matrix"][geneset_index[i],geneset_index[j]],RW_dict["matrix"][geneset_index[j],geneset_index[i]])

        # Obtain labels
        if show_labels:
            logging.info("show symbols")
            converter = Converter(mapping_table_file="../../primary_data/entrez_name.tsv")
            labels = converter.entrez2symbol(item)
        else:
            labels=[i for i in item]


        Z = linkage(X, 'single')
        fig = plt.figure(figsize=(25, 10))
        dn = dendrogram(Z,labels=labels)
        plt.savefig(outpath+setname+".pdf",f="pdf")

        fig =plt.figure(figsize=(25, 10))
        plt.imshow(1-X)
        g2=sns.heatmap(1-X,cmap="OrRd",square=True,xticklabels=labels,yticklabels=labels,cbar=True)
        g2.set_yticklabels(g2.get_yticklabels(), rotation = 0, fontsize = 6)
        g2.set_xticklabels(g2.get_xticklabels(), rotation = 90, fontsize = 6)

        plt.savefig(outpath+setname+"_matrix.pdf",f="pdf")

        subgraph=network.subgraph(item)
        layout = nx.spring_layout(subgraph)
        plt.figure(figsize=(13.5,5))
        plt.subplot(1,3,1)
        nx.draw_networkx(subgraph, nodelist=subgraph.nodes(), pos=layout, node_color='#636363',fontsize=4,  node_size = 30, edge_color='#bdbdbd' )
        #nx.draw_networkx_labels(subgraph, nodelist=subgraph.nodes(), pos=layout, fontsize=4, node_color='#636363', node_size = 30, edge_color='#bdbdbd' )
        plt.savefig(outpath+setname+"_graph.pdf", f="pdf")



        if show_results:
            painter_rw = paint.Painter_RW(
                network, output1.output, setname, RW_dict, item)
            logging.info("Painting matrix")
            painter_rw.plot_matrix()

    #output1.save_output_summary()

    if show_results:
        paint.paint_datasets_stats(
            output1.output_table, output1.output, "RWR")


def hierarchical_clustering_SP(network_file: "network file, use a network with weights",
                                geneset_file: "GMT geneset file",
                                distance_matrix_filename: "hdf5 obtained calculating a shortest path matrix",
                                outpath:"outpath",
                                setname: "Geneset to analyse" = None,
                                show_results: "image of the dendrogram" = False,
                                show_labels=False):

    """
        The script performs the hierarchical clustering of the RW matrix of a geneset
    """

    network = __load_network(network_file)
    network = nx.Graph(network.subgraph(
          max(nx.connected_components(network), key=len)))

    geneset = __load_geneset(geneset_file, setname)

    diz = {}
    diz["nodes"], diz["matrix"] = __read_distance_matrix(distance_matrix_filename)

    #create output
    setnames = [key for key in geneset.keys()]
    output1 = out.Output(network_file, outpath, prefix,
                         "hierarchical clustering", geneset_file, setnames)
    logging.info("Output-folder= " + output1.output)
    output1.set_diffusion_matrix(distance_matrix_filename)

    # Clustering for each dataset
    for setname, item in geneset.items():
        logging.info("\nSetname: " + setname)

        nodes_names = [i.decode("utf-8") for i in diz["nodes"]]
        item = set(item).intersection(set(nodes_names))

        #mapping each gene to the ones of the matrix
        geneset_index=[diz["nodes"].index(i.encode('utf-8')) for i in item if i.encode('utf-8') in diz["nodes"]]
        logging.info("There are %d nodes to be clustered" %len(geneset_index))

        # Allocate a matrix, it is going to be the "distance" matrix
        X = np.ones((len(geneset_index),len(geneset_index)))
        for i in range(len(geneset_index)):
            for j in range(len(geneset_index)):
                if i!=j:
                    # we use the distance between the two genes
                    X[i,j]=max(diz["matrix"][geneset_index[i],geneset_index[j]],diz["matrix"][geneset_index[j],geneset_index[i]])

        # Obtain labels
        if show_labels:
            logging.info("show symbols")
            converter = Converter(mapping_table_file="../../primary_data/entrez_name.tsv")
            labels = converter.entrez2symbol(item)
        else:
            labels=[i for i in item]


        Z = linkage(X, 'average', optimal_ordering=True)
        print(Z)


        max_d=20.5

        c=fcluster(Z,max_d,criterion="distance")
        print("cluster")

        clusters=set(c)
        print(c)
        c_number=0
        for clust in clusters:
            index_in_clust,=np.where(c == clust)
            nodes_in_clust = [diz["nodes"][i].decode('utf-8') for i in index_in_clust]
            output1.add_GMT_entry("Cluster_%d" %c_number, "cluster_SP", nodes_in_clust)
            c_number+=1

        output1.create_GMT_output()

        fig = plt.figure(figsize=(25, 10))
        fancy_dendrogram(
            Z,
            #truncate_mode='lastp',
            #p=12,
            leaf_rotation=90.,
            leaf_font_size=12.,
            #show_contracted=True,
            annotate_above=10,
            max_d=max_d,  # plot a horizontal cut-off line
        )
        #plt.show()


        #fig = plt.figure(figsize=(25, 10))
        #dn = dendrogram(Z,labels=labels)
        plt.savefig(output1.output+setname+".pdf",f="pdf")

        #fig =plt.figure(figsize=(25, 10))
        #plt.imshow(X)
        #g2=sns.heatmap(X,cmap="OrRd",square=True,xticklabels=labels,yticklabels=labels,cbar=True)
        #g2.set_yticklabels(g2.get_yticklabels(), rotation = 0, fontsize = 6)
        #g2.set_xticklabels(g2.get_xticklabels(), rotation = 90, fontsize = 6)

        #plt.savefig(output1.output_path+output1.output_folder+"/"+setname+"_matrix.pdf",f="pdf")

        #subgraph=network.subgraph(item)
        #layout = nx.spring_layout(subgraph)

        #plt.figure(figsize=(13.5,5))
        #plt.subplot(1,3,1)
        #nx.draw_networkx(subgraph, nodelist=subgraph.nodes(), pos=layout, node_color='#636363',fontsize=4,  node_size = 30, edge_color='#bdbdbd' )
        #nx.draw_networkx_labels(subgraph, nodelist=subgraph.nodes(), pos=layout, fontsize=4, node_color='#636363', node_size = 30, edge_color='#bdbdbd' )
        #plt.savefig(outpath+setname+"_graph.pdf", f="pdf")



        if show_results:
            painter_rw = paint.Painter_RW(
                network, output1.output, setname, diz, item)
            logging.info("Painting matrix")
            painter_rw.plot_matrix()

    #output1.save_output_summary()

    if show_results:
        paint.paint_datasets_stats(
            output1.output_table, output1.output, "SP")
