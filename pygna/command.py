"""PyGNA commands main module
"""
import logging
import numpy as np
import pygna.output as out
import networkx as nx
import pygna.elaborators as pe
import pygna.reading_class as rc
import pygna.statistical_test as st
import pygna.statistical_comparison as sc
import pygna.diagnostic as diagnostic
import pygna.painter as paint
import pygna.utils as utils
import pygna.statistical_diffusion as sd
import scipy
import scipy.linalg.interpolative
import itertools
import tables


def read_distance_matrix(distance_matrix_filename, in_memory=False):
    """
    Reads the large matrix. Uses a hdf5 file to work with it

    :param distance_matrix_filename: str, the file to read
    :param in_memory: bool, if the table must be kept in memory
    """
    nodes, data = rc.ReadDistanceMatrix(distance_matrix_filename, in_memory).get_data()
    return nodes, data


################################################################################
######### SUMMARY ##############################################################
################################################################################


def network_summary(network_file: "network file",
                    text_output: "output text file for the summary",
                    degree_figure_file: "pdf or png file for the degree distribution",
                    c_components_figure_file: "pdf or png file for the connected components distribution",
                    geneset_input_file: "geneset file" = None,
                    setname: "specify a single geneset" = None):
    """
    This function saves the principal info of a graph:
    - network properties
    - degree distribution
    - connected components diagnostic

    If a geneset/setname is passed to the function, the properties of
    the subgraph are evaluated
    """
    logging.info("Evaluating network summary, please wait")
    network = rc.ReadTsv(network_file).get_network()

    if geneset_input_file:
        if not setname:
            logging.error('Missing setname name, specify  a unique setname')
        else:
            geneset = rc.ReadGmt(geneset_input_file).get_geneset(setname)
            for setname, item in geneset.items():
                graph = nx.subgraph(network, item)
                out.write_graph_summary(graph, text_output, setname + " on " + network_file)
                diagnostic.plot_connected_components(nx.connected_components(graph), c_components_figure_file)
                diagnostic.plot_degree(nx.degree(graph), degree_figure_file)

    else:
        out.write_graph_summary(network, text_output, network_file)
        diagnostic.plot_connected_components(nx.connected_components(network), c_components_figure_file)
        diagnostic.plot_degree(nx.degree(network), degree_figure_file)
    logging.info("Network summary completed")


################################################################################
######### SINGLE SET ANALYSES ##################################################
################################################################################


def test_topology_total_degree(
    network_file: "network file",
    geneset_file: "GMT geneset file",
    output_table: "output results table, use .csv extension",
    setname: "Geneset to analyse" = None,
    size_cut: "removes all genesets with a mapped length < size_cut" = 20,
    number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
    cores: "Number of cores for the multiprocessing" = 1,
    results_figure: "barplot of results, use pdf or png extension" = None,
    diagnostic_null_folder: "plot null distribution, pass the folder where all the figures are going to be saved "
                            "(one for each dataset)" = None):
    """
        Performs the analysis of total degree of the .

        It computes a p-value for the ratio of total degree
        of the geneset being bigger than the one expected by chance
        for a geneset of the same size.
    """
    logging.info("Evaluating the test topology total degree, please wait")
    network = rc.ReadTsv(network_file).get_network()

    geneset = rc.ReadGmt(geneset_file).get_geneset(setname)

    setnames = [key for key in geneset.keys()]

    # Generate output
    output1 = out.Output(network_file, output_table, "topology_total_degree", geneset_file, setnames)
    logging.info("Results file = " + output1.output_table_results)

    # Create table
    output1.create_st_table_empirical()
    st_test = st.StatisticalTest(st.geneset_total_degree_statistic, network)

    for setname, item in geneset.items():
        # Geneset smaller than size cut are not taken into consideration
        if len(item) > size_cut:
            item = set(item)
            observed, pvalue, null_d, n_mapped, n_geneset = st_test.empirical_pvalue(item,
                                                                                     max_iter=number_of_permutations,
                                                                                     alternative="greater",
                                                                                     cores=cores)
            logging.info("Setname:" + setname)
            if n_mapped < size_cut:
                logging.info("%s removed from results since nodes mapped are < %d" % (setname, size_cut))
            else:
                logging.info("Observed: %g p-value: %g" % (observed, pvalue))
                # TODO Check line below
                logging.info("Null mean: %g null variance: %g".format(np.mean(null_d), np.var(null_d)))
                output1.update_st_table_empirical(setname, n_mapped, n_geneset, number_of_permutations, observed,
                                                  pvalue, np.mean(null_d), np.var(null_d))
                if diagnostic_null_folder:
                    diagnostic.plot_null_distribution(null_d, observed, diagnostic_null_folder + setname +
                                                      '_total_degree_null_distribution.pdf', setname=setname)
    output1.close_temporary_table()
    if results_figure:
        paint.paint_datasets_stats(output1.output_table_results, results_figure, alternative='greater')
    logging.info("Test topology total degree completed")


def test_topology_internal_degree(
    network_file: "network file",
    geneset_file: "GMT geneset file",
    output_table: "output results table, use .csv extension",
    setname: "Geneset to analyse" = None,
    size_cut: "removes all genesets with a mapped length < size_cut" = 20,
    number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
    cores: "Number of cores for the multiprocessing" = 1,
    results_figure: "barplot of results, use pdf or png extension" = None,
    diagnostic_null_folder: "plot null distribution, pass the folder where all the figures are going to be saved "
                            "(one for each dataset)" = None,
):
    """
        Performs the analysis of internal degree.
        It computes a p-value for the ratio of internal degree
        of the geneset being bigger than the one expected by chance
        for a geneset of the same size.
    """

    network = rc.ReadTsv(network_file).get_network()
    geneset = rc.ReadGmt(geneset_file).get_geneset(setname)
    setnames = [key for key in geneset.keys()]
    output1 = out.Output(network_file, output_table, "topology_internal_degree", geneset_file, setnames)
    logging.info("Results file = " + output1.output_table_results)
    output1.create_st_table_empirical()
    st_test = st.StatisticalTest(st.geneset_internal_degree_statistic, network)
    for setname, item in geneset.items():
        item = set(item)
        if len(item) > size_cut:
            observed, pvalue, null_d, n_mapped, n_geneset = st_test.empirical_pvalue(item,
                                                                                     max_iter=number_of_permutations,
                                                                                     alternative="greater", cores=cores)

            logging.info("Setname:" + setname)
            if n_mapped < size_cut:
                logging.info("%s remove from results since nodes mapped are < %d" % (setname, size_cut))
            else:
                logging.info("Observed: %g p-value: %g" % (observed, pvalue))
                output1.update_st_table_empirical(setname, n_mapped, n_geneset, number_of_permutations, observed,
                                                  pvalue, np.mean(null_d), np.var(null_d))

                if diagnostic_null_folder:
                    diagnostic.plot_null_distribution(null_d, observed, diagnostic_null_folder + setname +
                                                      '_internal_degree_null_distribution.pdf', setname=setname)
    output1.close_temporary_table()
    if results_figure:
        paint.paint_datasets_stats(output1.output_table_results, results_figure, alternative='greater')


def test_topology_rwr(
    network_file: "network file, use a network with weights",
    geneset_file: "GMT geneset file",
    rwr_matrix_filename: "hdf5 RWR matrix obtained with pygna ",
    output_table: "output results table, use .csv extension",
    setname: "Geneset to analyse" = None,
    size_cut: "removes all genesets with a mapped length < size_cut" = 20,
    number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
    cores: "Number of cores for the multiprocessing" = 1,
    in_memory: "set if you want the large matrix to be read in memory" = False,
    results_figure: "barplot of results, use pdf or png extension" = None,
    diagnostic_null_folder: "plot null distribution, pass the folder where all the figures are going to be saved "
                            "(one for each dataset)" = None,
):
    """
        Performs the analysis of random walk probabilities.
        Given the RWR matrix ,
        it compares the probability of walking between the genes in the geneset
        compared to those of walking between the nodes
        of a geneset with the same size
    """

    network = rc.ReadTsv(network_file).get_network()
    network = nx.Graph(network.subgraph(max(nx.connected_components(network), key=len)))
    geneset = rc.ReadGmt(geneset_file).get_geneset(setname)
    rw_dict = {"nodes": read_distance_matrix(rwr_matrix_filename, in_memory=in_memory)[0],
               "matrix": read_distance_matrix(rwr_matrix_filename, in_memory=in_memory)[1]}

    setnames = [key for key in geneset.keys()]
    output1 = out.Output(network_file, output_table, "topology_rwr", geneset_file, setnames)

    logging.info("Results file = " + output1.output_table_results)
    output1.create_st_table_empirical()
    st_test = st.StatisticalTest(st.geneset_RW_statistic, network, rw_dict)

    for setname, item in geneset.items():
        item = set(item)
        if len(item) > size_cut:
            # test
            observed, pvalue, null_d, n_mapped, n_geneset = st_test.empirical_pvalue(item,
                                                                                     max_iter=number_of_permutations,
                                                                                     alternative="greater", cores=cores)
            logging.info("Setname:" + setname)
            if n_mapped < size_cut:
                logging.info("%s remove from results since nodes mapped are < %d" % (setname, size_cut))
            else:
                logging.info("Observed: %g p-value: %g" % (observed, pvalue))
                if diagnostic_null_folder:
                    diagnostic.plot_null_distribution(null_d, observed, diagnostic_null_folder + setname +
                                                      '_rwr_null_distribution.pdf', setname=setname)
                # saving output
                output1.update_st_table_empirical(setname, n_mapped, n_geneset, number_of_permutations, observed,
                                                  pvalue, np.mean(null_d), np.var(null_d))
        else:
            logging.info("%s removed from results since nodes mapped are < %d" % (setname, size_cut))

    output1.close_temporary_table()
    if results_figure:
        paint.paint_datasets_stats(output1.output_table_results, results_figure, alternative='greater')


def test_topology_module(
    network_file: "network file",
    geneset_file: "GMT geneset file",
    output_table: "output results table, use .csv extension",
    setname: "Geneset to analyse" = None,
    size_cut: "removes all genesets with a mapped length < size_cut" = 20,
    number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
    cores: "Number of cores for the multiprocessing" = 1,
    output_lcc: "for creating a GMT file with the LCC lists pass a gmt filename" = None,
    results_figure: "barplot of results, use pdf or png extension" = None,
    diagnostic_null_folder: "plot null distribution, pass the folder where all the figures are going to be saved "
                            "(one for each dataset)" = None,
):
    """
        Performs geneset network topology module analysis.

        It computes a p-value for the largest connected component
        of the geneset being bigger than the one expected by chance
        for a geneset of the same size.
    """
    network = rc.ReadTsv(network_file).get_network()
    geneset = rc.ReadGmt(geneset_file).get_geneset(setname)

    setnames = [key for key in geneset.keys()]
    output1 = out.Output(network_file, output_table, "topology_module", geneset_file, setnames)
    logging.info("Results file = " + output1.output_table_results)
    output1.create_st_table_empirical()

    st_test = st.StatisticalTest(st.geneset_module_statistic, network)
    for setname, item in geneset.items():
        item = set(item)
        if len(item) > size_cut:
            if output_lcc:
                module = nx.subgraph(network, item)
                if len(module.nodes) > 0:
                    lcc = sorted(list(nx.connected_components(module)), key=len, reverse=True)[0]
                else:
                    lcc = []
                output1.add_GMT_entry(setname, "topology_module", lcc)

            observed, pvalue, null_d, n_mapped, n_geneset = st_test.empirical_pvalue(item,
                                                                                     max_iter=number_of_permutations,
                                                                                     alternative="greater", cores=cores)
            logging.info("Setname:" + setname)
            if n_mapped < size_cut:
                logging.info("%s remove from results since nodes mapped are < %d" % (setname, size_cut))
            else:
                logging.info("Observed: %g p-value: %g" % (observed, pvalue))
                output1.update_st_table_empirical(setname, n_mapped, n_geneset, number_of_permutations, observed,
                                                  pvalue, np.mean(null_d), np.var(null_d))
                if diagnostic_null_folder:
                    diagnostic.plot_null_distribution(null_d, observed, diagnostic_null_folder + setname +
                                                      '_module_null_distribution.pdf', setname=setname)
    output1.close_temporary_table()
    if output_lcc:
        output1.create_GMT_output(output_lcc)

    if results_figure:
        paint.paint_datasets_stats(output1.output_table_results, results_figure, alternative='greater')


def test_topology_sp(
    network_file: "network file",
    geneset_file: "GMT geneset file",
    distance_matrix_filename: "distance hdf5 matrix file generated by pygna",
    output_table: "output results table, use .csv extension",
    setname: "Geneset to analyse" = None,
    size_cut: "removes all genesets with a mapped length < size_cut" = 20,
    number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
    cores: "Number of cores for the multiprocessing" = 1,
    in_memory: "set if you want the large matrix to be read in memory" = False,
    results_figure: "barplot of results, use pdf or png extension" = None,
    diagnostic_null_folder: "plot null distribution, pass the folder where all the figures are going to be saved "
                            "(one for each dataset)" = None,
):
    """
        Performs geneset network topology shortest path analysis.

        It computes a p-value for the average shortest path length
        of the geneset being smaller than expected by chance
        for a geneset of the same size.

    """

    network = rc.ReadTsv(network_file).get_network()
    network = nx.Graph(network.subgraph(max(nx.connected_components(network), key=len)))

    geneset = rc.ReadGmt(geneset_file).get_geneset(setname)

    diz = {"nodes": read_distance_matrix(distance_matrix_filename, in_memory=in_memory)[0],
           "matrix": read_distance_matrix(distance_matrix_filename, in_memory=in_memory)[1]}
    diz["matrix"] = diz["matrix"] + np.transpose(diz["matrix"])
    np.fill_diagonal(diz["matrix"], float("inf"))
    setnames = [key for key in geneset.keys()]

    output1 = out.Output(network_file, output_table, "topology_sp", geneset_file, setnames)
    logging.info("Results file = " + output1.output_table_results)
    output1.create_st_table_empirical()
    st_test = st.StatisticalTest(st.geneset_localisation_statistic, network, diz)

    for setname, item in geneset.items():

        item = set(item)
        if len(item) > size_cut:
            logging.info("Setname:" + setname)
            observed, pvalue, null_d, n_mapped, n_geneset = st_test.empirical_pvalue(item, cores=cores,
                                                                                     max_iter=number_of_permutations)
            logging.info("Observed: %g p-value: %g" % (observed, pvalue))

            output1.update_st_table_empirical(setname, n_mapped, n_geneset, number_of_permutations, observed, pvalue,
                                              np.mean(null_d), np.var(null_d))
            if diagnostic_null_folder:
                diagnostic.plot_null_distribution(null_d, observed, diagnostic_null_folder + setname +
                                                  '_sp_null_distribution.pdf', setname=setname, alternative="less")
        else:
            logging.info("%s remove from results since nodes mapped are < %d" % (setname, size_cut))
    output1.close_temporary_table()
    if results_figure:
        paint.paint_datasets_stats(output1.output_table_results, results_figure, alternative='less')


################################################################################
######### Diffusion test #######################################################
################################################################################
def test_diffusion_hotnet(network_file: "network file, use a network with weights",
                          geneset_file: "csv geneset file",
                          rwr_matrix_filename: "hdf5 RWR matrix obtained with pygna ",
                          output_table: "output results table, use .csv extension",
                          name_column: "Column to use as name (default is deseq2)" = "gene_name",
                          weight_column: "Column to use as weight (default is deseq2)" = "stat",
                          filter_column: "Column used to define the significant genes (default is deseq2)" = "padj",
                          filter_condition: "Condition for significance" = "less",
                          filter_threshold: "threshold for significance" = 0.01,
                          normalise: 'pass this flag for using only positive values in the analysis' = False,
                          size_cut: "removes all genesets with a mapped length < size_cut" = 20,
                          number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
                          cores: "Number of cores for the multiprocessing" = 1,
                          in_memory: "set if you want the large matrix to be read in memory" = False,
                          ):
    """
        Performs the analysis of random walk applying the weights of an upstream analysis.
        Given a csv file the user needs to specify the columns of interest and
        the threshold of significance.
        For the analysis the StatisticalDiffusion is used with hotnet_diffusion_statistic
        function.
    """

    # Reading network file
    network = rc.ReadTsv(network_file).get_network()
    network = nx.Graph(network.subgraph(max(nx.connected_components(network), key=len)))

    # Read geneset
    table = rc.ReadCsv(geneset_file, column_to_fill=name_column).get_data()
    if len(table.columns) < 2:
        logging.error("Error: the function takes a csv file as input, the read file has less than 2 columns, "
                      "check that the table is comma separated")

    # Filter table for significant genes
    table[name_column] = table[name_column].fillna(0).apply(str)
    table = pe.TableElaboration.clean_table(table=table, stat_col=weight_column)
    geneset = utils.filter_table(table, filter_column=filter_column, alternative=filter_condition,
                                 threshold=filter_threshold)[name_column]
    if normalise:
        table[weight_column] = np.abs(table[weight_column].values)

    if len(geneset) < size_cut:
        logging.error('The number of significant genes is lower than %d. \
                    \n Change size_cut if necessary' % size_cut)

    # Read RWR matrix
    rw_dict = {"nodes": read_distance_matrix(rwr_matrix_filename, in_memory=in_memory)[0],
               "matrix": read_distance_matrix(rwr_matrix_filename, in_memory=in_memory)[1]}
    # setting output
    output1 = out.Output(network_file, output_table, "diffusion", geneset_file, geneset_file)
    output1.create_st_table_empirical()
    logging.info("Results file = " + output1.output_table_results)

    # initialising test
    st_test = sd.DiffusionTest(sd.hotnet_diffusion_statistic, rw_dict["nodes"], rw_dict["matrix"], table,
                               names_col=name_column, weights_col=weight_column)

    observed, pvalue, null_d, n_mapped, n_geneset = st_test.empirical_pvalue(geneset, max_iter=number_of_permutations,
                                                                             alternative="greater", cores=cores)
    if n_mapped < size_cut:
        logging.info("Results removed, since nodes mapped are < %d" % size_cut)
    else:
        logging.info("Observed: %g p-value: %g" % (observed, pvalue))
        output1.update_st_table_empirical(geneset_file, n_mapped, n_geneset, number_of_permutations, observed, pvalue,
                                          np.mean(null_d), np.var(null_d))
    output1.close_temporary_table()

    # if results_figure:
    #    paint.paint_diffusion_matrix(output1.output_table_results, results_figure, alternative='greater',)


################################################################################
######### ASSOCIATIONS and COMPARISONS #########################################
################################################################################

def test_association_sp(
    network_file: "network file",
    file_geneset_a: "GMT geneset file, if it's the only parameter passed the analysis is gonna be run on all the "
                    "couples of datasets, otherwise specify the other files and setnames",
    distance_matrix_filename: "distance matrix file generated by pygna",
    output_table: "output results table, use .csv extension",
    setname_a: "Geneset A to analyse" = None,
    file_geneset_b: "GMT geneset file" = None,
    setname_b: "Geneset B to analyse" = None,
    size_cut: "removes all genesets with a mapped length < size_cut" = 20,
    keep: "if true, keeps the geneset B not permuted" = False,
    cores: "Number of cores for the multiprocessing" = 1,
    in_memory: "set if you want the large matrix to be read in memory" = False,
    number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
    results_figure: "barplot of results, use pdf or png extension" = None,
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

    network = rc.ReadTsv(network_file).get_network()
    network = nx.Graph(network.subgraph(max(nx.connected_components(network), key=len)))

    # Read matrix
    sp_diz = {"nodes": read_distance_matrix(distance_matrix_filename, in_memory=in_memory)[0],
              "matrix": read_distance_matrix(distance_matrix_filename, in_memory=in_memory)[1]}
    # TODO check why "in memory"= True results are different from "in memory"= False
    sp_diz["matrix"] = sp_diz["matrix"] + np.transpose(sp_diz["matrix"])
    np.fill_diagonal(sp_diz["matrix"], np.inf)

    # Managing the different genesets
    if setname_a and setname_b is None and file_geneset_b is None:
        logging.error(" this analysis requires at least two genesets ")

    geneset_a = rc.ReadGmt(file_geneset_a).get_geneset(setname_a)
    if file_geneset_b:
        geneset_b = rc.ReadGmt(file_geneset_b).get_geneset(setname_b)
    else:
        if setname_b:
            geneset_b = rc.ReadGmt(file_geneset_a).get_geneset(setname_b)
        else:
            geneset_b = None

    st_comparison = sc.StatisticalComparison(sc.comparison_shortest_path, network, diz=sp_diz, n_proc=cores)

    if not geneset_b:  # Analysis of genesets inside a single file
        logging.info("Analysing all the sets in " + file_geneset_a)
        setnames = [key for key in geneset_a.keys()]

        # Creating the output table
        output1 = out.Output(network_file, output_table, analysis_name_str, file_geneset_a, setnames)
        logging.info("Results file = " + output1.output_table_results)
        output1.create_comparison_table_empirical()

        for pair in itertools.combinations(setnames, 2):
            if len(set(geneset_a[pair[0]])) > size_cut and len(set(geneset_a[pair[1]])) > size_cut:
                logging.info("Analysing " + str(pair[0]) + " and " + str(pair[1]))

                n_overlaps = len(set(geneset_a[pair[0]]).intersection(set(geneset_a[pair[1]])))
                observed, pvalue, null_d, a_mapped, b_mapped = st_comparison.comparison_empirical_pvalue(
                    set(geneset_a[pair[0]]), set(geneset_a[pair[1]]), max_iter=number_of_permutations, keep=keep)
                # Save the results

                output1.update_comparison_table_empirical(pair[0], pair[1], len(set(geneset_a[pair[0]])), a_mapped,
                                                          len(set(geneset_a[pair[1]])), b_mapped, n_overlaps,
                                                          number_of_permutations, observed, pvalue, np.mean(null_d),
                                                          np.var(null_d))
            else:
                logging.warning("Geneset A has %d terms and Geneset B has %d terms. \
                \nOne of them is too short, analysis not done" % (len(set(geneset_a[pair[0]])),
                                                                  len(set(geneset_a[pair[1]]))))

    else:  # Analysis of genesets into two different gmt files

        logging.info("geneset_a contains %d sets", (len(geneset_a)))
        sets_a = [key for key in geneset_a.keys()]
        logging.info("geneset_b contains %d sets", (len(geneset_b)))
        sets_b = [key for key in geneset_b.keys()]
        output1 = out.Output(network_file, output_table, analysis_name_str, file_geneset_a, sets_a, file_geneset_b,
                             sets_b)
        logging.info("Results file = " + output1.output_table_results)
        output1.create_comparison_table_empirical()
        for set_A, item_A in geneset_a.items():
            for set_B, item_B in geneset_b.items():
                n_overlaps = len(set(item_A).intersection(set(item_B)))
                if len(item_A) > size_cut and len(item_B) > size_cut:
                    observed, pvalue, null_d, a_mapped, b_mapped = st_comparison.comparison_empirical_pvalue(
                        set(item_A), set(item_B), max_iter=number_of_permutations, keep=keep)

                    logging.info("Observed: %g p-value: %g" % (observed, pvalue))
                    output1.update_comparison_table_empirical(set_A, set_B, len(set(item_A)), a_mapped,
                                                              len(set(item_B)), b_mapped, n_overlaps,
                                                              number_of_permutations, observed, pvalue, np.mean(null_d),
                                                              np.var(null_d))
    output1.close_temporary_table()
    if results_figure:
        paint.paint_comparison_matrix(output1.output_table_results, results_figure)


def test_association_rwr(
    network_file: "network file",
    file_geneset_a: "GMT geneset file",
    rwr_matrix_filename: ".hdf5 file with the RWR matrix obtained by pygna",
    output_table: "output results table, use .csv extension",
    setname_a: "Geneset A to analyse" = None,
    file_geneset_b: "GMT geneset file" = None,
    setname_b: "Geneset B to analyse" = None,
    size_cut: "removes all genesets with a mapped length < size_cut" = 20,
    keep: "if true, keeps the geneset B unpermuted" = False,
    cores: "Number of cores for the multiprocessing" = 1,
    in_memory: "set if you want the large matrix to be read in memory" = False,
    number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
    results_figure: "heatmap of results" = None,
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

    network = rc.ReadTsv(network_file).get_network()
    network = nx.Graph(network.subgraph(max(nx.connected_components(network), key=len)))

    # Read datasets
    if setname_a and setname_b is None and file_geneset_b is None:
        logging.error(" this analysis requires at least two genesets ")

    rw_dict = {"nodes": read_distance_matrix(rwr_matrix_filename, in_memory=in_memory)[0],
               "matrix": read_distance_matrix(rwr_matrix_filename, in_memory=in_memory)[1]}

    geneset_a = rc.ReadGmt(file_geneset_a).get_geneset(setname_a)
    if file_geneset_b:
        geneset_b = rc.ReadGmt(file_geneset_b).get_geneset(setname_b)
    else:
        if setname_b:
            geneset_b = rc.ReadGmt(file_geneset_a).get_geneset(setname_b)
        else:
            geneset_b = None

    st_comparison = sc.StatisticalComparison(sc.comparison_random_walk, network, n_proc=cores, diz=rw_dict)
    if not geneset_b:
        logging.info("Analysing all the sets in " + file_geneset_a)
        setnames = [key for key in geneset_a.keys()]
        output1 = out.Output(network_file, output_table, analysis_name_str, file_geneset_a, setnames)
        logging.info("Results file = " + output1.output_table_results)
        output1.create_comparison_table_empirical()

        for pair in itertools.combinations(setnames, 2):
            if len(set(geneset_a[pair[0]])) > size_cut and len(set(geneset_a[pair[1]])) > size_cut:
                logging.info("Analysing " + str(pair[0]) + " and " + str(pair[1]))
                n_overlaps = len(set(geneset_a[pair[0]]).intersection(set(geneset_a[pair[1]])))
                observed, pvalue, null_d, a_mapped, b_mapped = st_comparison.comparison_empirical_pvalue(
                    set(geneset_a[pair[0]]), set(geneset_a[pair[1]]), max_iter=number_of_permutations,
                    alternative="greater", keep=keep)

                output1.update_comparison_table_empirical(pair[0], pair[1], len(set(geneset_a[pair[0]])), a_mapped,
                                                          len(set(geneset_a[pair[1]])), b_mapped, n_overlaps,
                                                          number_of_permutations, observed, pvalue, np.mean(null_d),
                                                          np.var(null_d))
            else:
                logging.warning("Geneset A has %d terms and Geneset B has %d terms. \
                \nOne of them is too short, analysis not done" % (len(set(geneset_a[pair[0]])),
                                                                  len(set(geneset_a[pair[1]]))))

    else:
        logging.info("geneset_a contains %d sets" % (len(geneset_a)))
        sets_a = [key for key in geneset_a.keys()]
        logging.info("Setnames in A: " + str(sets_a))
        logging.info("geneset_b contains %d sets" % (len(geneset_b)))
        sets_b = [key for key in geneset_b.keys()]
        logging.info("Setnames in B: " + str(sets_b))
        output1 = out.Output(network_file, output_table, analysis_name_str, file_geneset_a, sets_a, file_geneset_b,
                             sets_b)
        logging.info("Results file = " + output1.output_table_results)
        output1.create_comparison_table_empirical()

        for set_A, item_A in geneset_a.items():
            for set_B, item_B in geneset_b.items():

                if len(item_A) > size_cut and len(item_B) > size_cut:
                    logging.info("Analysing " + str(set_A) + " and " + str(set_B))
                    n_overlaps = len(set(item_A).intersection(set(item_B)))
                    observed, pvalue, null_d, a_mapped, b_mapped = st_comparison.comparison_empirical_pvalue(
                        set(item_A), set(item_B), max_iter=number_of_permutations, alternative="greater", keep=keep)
                    logging.info("Observed: %g p-value: %g" % (observed, pvalue))

                    output1.update_comparison_table_empirical(set_A, set_B, len(set(item_A)), a_mapped,
                                                              len(set(item_B)), b_mapped, n_overlaps,
                                                              number_of_permutations, observed, pvalue, np.mean(null_d),
                                                              np.var(null_d))
                else:
                    logging.warning("Geneset A has %d terms and Geneset B has %d terms. \
                    \nOne of them is too short, analysis not done" % (len(set(item_A)), len(set(item_B))))

    output1.close_temporary_table()
    if results_figure:
        paint.paint_comparison_matrix(output1.output_table_results, results_figure, rwr=True)


################################################################################
######### BUILDING FUNCTIONS ###################################################
################################################################################

def build_distance_matrix(
    network_file: "network file",
    output_file: "distance matrix output file, use .hdf5",
    giant_component_only: "compute the shortest paths only for nodes in the giant component" = True,
):
    """
        Build a shortest path distance matrix for a given network.
        Matrix can be saved as a lm.txt file or a .hdf5 one.
    """
    logging.info("Converting distance matrix, please wait...")
    network = rc.ReadTsv(network_file).get_network()
    if giant_component_only:
        network = network.subgraph(max(nx.connected_components(network), key=len))

    distance_matrix = nx.all_pairs_shortest_path_length(network)

    nodes = list(network.nodes())

    if output_file.endswith(".hdf5"):
        with tables.open_file(output_file, mode="w") as hdf5_file:

            # create a hdf5 file with two objects:
            # - one is the nodes array,
            hdf5_nodes = hdf5_file.create_array(hdf5_file.root, "nodes", nodes)
            # -  the other is the shortest path distance matrix
            hdf5_data = hdf5_file.create_array(hdf5_file.root, "matrix", np.zeros((len(nodes), len(nodes))))
            for node_i, k in distance_matrix:
                i = nodes.index(node_i)
                for node_j, sp in k.items():
                    j = nodes.index(node_j)
                    if j >= i:  # saves only the upper triangular matrix
                        hdf5_data[i, j] = sp
        hdf5_file.close()
    else:
        logging.error("Pass an hd5f file")


def build_rwr_diffusion(
    network_file: "network file",
    beta=0.85,
    output_file: "distance matrix output file (use .hdf5) " = None,
):
    """
        Build the RWR_diffusion_matrix
    """

    network = rc.ReadTsv(network_file).get_network()
    network = network.subgraph(max(nx.connected_components(network), key=len))
    nodes = list(network.nodes())

    logging.info("Beginning to calculate RWR matrix")

    a = nx.adjacency_matrix(network)
    k = 1 / np.array(list(dict(network.degree()).values()))
    d = scipy.sparse.dia_matrix((k, [0]), shape=a.shape)
    a = a.dot(d)
    n = np.shape(a)[1]

    if output_file.endswith(".hdf5"):
        with tables.open_file(output_file, mode="w") as hdf5_file:
            # create a hdf5 file with two objects:
            # - one is the nodes array,
            hdf5_file.create_array(hdf5_file.root, "nodes", nodes)
            # -  the other is the RWR matrix
            hdf5_file.create_array(hdf5_file.root, "matrix", beta * np.linalg.inv(np.eye(n) - (1.0 - beta) * a))
            logging.info("Saving network")
            hdf5_file.close()
    else:
        return beta * np.linalg.inv(np.eye(n) - (1.0 - beta) * a)


def network_graphml(
    network_file: "network file",
    geneset_file: "geneset file",
    output_file: "graphml file for network for visualisation",
    setname: "set name" = None,
    giant_component_only: "saves only the giant component of the network" = True,
    minimal: 'saves only the minimal graph' = False,
):
    """
    This function generates a graphml file with nodes annotation.
    Given a geneset, with k setnames, each node has k False/True
    annotations for each set.

    Warning: without minimal, this function saves the full network.
    The minimal graph saves only the nodes in the geneset and those that
    connect them with a shortest path.
    """

    network = rc.ReadTsv(network_file).get_network()
    geneset = rc.ReadGmt(geneset_file).get_geneset(setname)

    if giant_component_only:
        network = network.subgraph(max(nx.connected_components(network), key=len))

    dict_nodes = {}
    if minimal:
        for setname in geneset:
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
            for n in new_network_minimal.nodes():
                if n in geneset[setname]:
                    dict_nodes[n] = True
                else:
                    dict_nodes[n] = False
            nx.set_node_attributes(new_network_minimal, dict_nodes, setname)

    else:
        for setname in geneset:
            for n in network.nodes():
                if n in geneset[setname]:
                    dict_nodes[n] = True
                else:
                    dict_nodes[n] = False
            nx.set_node_attributes(network, dict_nodes, setname)

    nx.write_graphml(network, output_file)
