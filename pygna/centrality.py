import argh
import logging
import networkx as nx
import pygna.reading_class as rc
import pygna.output as out
import pygna.statistical_test as st
import pygna.painter as paint
import pygna.diagnostic as diagnostic
import numpy as np


def calculate_centrality(graph: nx.Graph, geneset: set, matrix: dict, observed_flag: bool = False) -> np.ndarray:
    """
    This function calculate the graph closeness centrality.
    It considers the whole graph and calculate the shortest path, then for each node in the graph calculate the node centrality as follows:

    :math:`node centrality = len(sp) -1 / tot_{sp}`

    where sp is the distance of the node with each other node and tot_sp is the total shortest paths for the whole graph.

    :param graph: The network to analyse
    :param geneset: the geneset to analyse
    :param matrix: The dictionary containing nodes and distance matrix
    """

    graph = nx.subgraph(graph, geneset)

    graph_centrality = list()
    for n in graph.nodes:
        matrix_id = matrix["nodes"].index(n)
        sp = matrix["matrix"][matrix_id]
        tot_sp = sum(sp)
        if tot_sp > 0:
            # Remove 1 because we are considering one node of the graph
            graph_centrality[n] = (len(sp) - 1) / tot_sp

    return np.asarray(graph_centrality)


def test_topology_centrality(
    network_file: "network file",
    geneset_file: "GMT geneset file",
    matrix: "The matrix with the SP for each node",
    output_table: "output results table, use .csv extension",
    setname: "Geneset to analyse" = None,
    size_cut: "removes all genesets with a mapped length < size_cut" = 20,
    number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
    cores: "Number of cores for the multiprocessing" = 1,
    results_figure: "barplot of results, use pdf or png extension" = None,
    diagnostic_null_folder: "plot null distribution, pass the folder where all the figures are going to be saved "
                            "(one for each dataset)" = None):

    logging.info("Evaluating the test topology total degree, please wait")
    network = rc.ReadTsv(network_file).get_network()
    geneset = rc.ReadGmt(geneset_file).get_geneset(setname)
    setnames = [key for key in geneset.keys()]
    # Generate output
    output1 = out.Output(network_file, output_table, "topology_total_degree", geneset_file, setnames)
    logging.info("Results file = " + output1.output_table_results)
    # Create table
    output1.create_st_table_empirical()
    st_test = st.StatisticalTest(calculate_centrality, network, matrix)
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
                logging.info("Null mean: %g null variance: %g".format(np.mean(null_d), np.var(null_d)))
                output1.update_st_table_empirical(setname, n_mapped, n_geneset, number_of_permutations, observed,
                                                  pvalue, np.mean(null_d), np.var(null_d))
                if diagnostic_null_folder:
                    diagnostic.plot_null_distribution(null_d, observed, diagnostic_null_folder + setname +
                                                      '_total_degree_null_distribution.pdf', setname=setname)
    output1.close_temporary_table()
    if results_figure:
        paint.paint_datasets_stats(output1.output_table_results, results_figure, alternative='greater')
    logging.info("Test topology CENTRALITY completed")


def main():
    """
    argh dispatch
    """
    argh.dispatch_commands([test_topology_centrality])


if __name__ == "__main__":
    """
    MAIN
    """
    main()
