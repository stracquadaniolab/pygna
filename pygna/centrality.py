import argh
import logging
import networkx as nx
import pygna.reading_class as rc
import pygna.output as out
import pygna.statistical_test as st
import pygna.painter as paint
import pygna.diagnostic as diagnostic
import pygna.command as cmd
import numpy as np


def average_closeness_centrality(graph: nx.Graph, geneset: set, diz: dict, observed_flag: bool = False) -> float:
    """
    This function calculates the average closeness centrality of a geneset. For
    a single node, the closeness centrality is defined as the inverse of the
    shortest path distance of the node from all the other nodes.

    Given a network with N nodes and a distance shortest path function between
    two nodes d(u,v) closeness centrality (u)= (N -1) / sum (v != u) d(u,v)

    :param graph: The network to analyse
    :param geneset: the geneset to analyse
    :param diz: The dictionary contains the nodes and a vector with the sum of shortest path for the whole network.
    """

    graph_centrality = []

    ids = [diz["nodes"].index(n) for n in geneset]
    graph_centrality = [(len(diz["nodes"]) - 1) / diz['vector'][idx] for idx in ids]

    return np.mean(graph_centrality)



def test_topology_centrality(
    network_file: "network file",
    geneset_file: "GMT geneset file",
    distance_matrix_filename: "The matrix with the SP for each node",
    output_table: "output results table, use .csv extension",
    setname: "Geneset to analyse" = None,
    size_cut: "removes all genesets with a mapped length < size_cut" = 20,
    number_of_permutations: "number of permutations for computing the empirical pvalue" = 500,
    cores: "Number of cores for the multiprocessing" = 1,
    in_memory: 'load hdf5 data onto memory' = False,):

    """
    This function calculates the average closeness centrality of a geneset.
    For a single node, the closeness centrality is defined as the inverse
    of the shortest path distance of the node from all the other nodes.

    """

    logging.info("Evaluating the test topology total degree, please wait")
    network = rc.ReadTsv(network_file).get_network()
    network = nx.Graph(network.subgraph(max(nx.connected_components(network), key=len)))
    geneset = rc.ReadGmt(geneset_file).get_geneset(setname)
    setnames = [key for key in geneset.keys()]


    diz = {"nodes": cmd.read_distance_matrix(distance_matrix_filename, in_memory=in_memory)[0],
           "matrix": cmd.read_distance_matrix(distance_matrix_filename, in_memory=in_memory)[1]}
    diz["matrix"] = diz["matrix"] + np.transpose(diz["matrix"])

    np.fill_diagonal(diz["matrix"], float(0))

    diz['vector'] = np.sum(diz["matrix"],axis = 0)

    # Generate output
    output1 = out.Output(network_file, output_table, "topology_centrality", geneset_file, setnames)
    logging.info("Results file = " + output1.output_table_results)
    # Create table
    output1.create_st_table_empirical()
    st_test = st.StatisticalTest(average_closeness_centrality, network, diz)

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
                output1.update_st_table_empirical(setname, n_mapped, n_geneset, number_of_permutations, observed,
                                                  pvalue, np.mean(null_d), np.var(null_d))

    output1.close_temporary_table()


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
