import networkx as nx


def calculate_centrality(graph: nx.Graph, matrix: dict) -> dict:
    """
    This function calculate the graph centrality.
    It considers the whole graph and calculate the shortest path, then for each node in the graph calculate the node centrality as follows:

    :math:`node centrality = len(sp) -1 / tot_{sp}`

    where sp is the distance of the node with each other node and tot_sp is the total shortest paths for the whole graph.

    :param graph: The network to analyse
    :param matrix: The dictionary containing nodes and distance matrix

    Example
    _______

    .. code-block:: python

            import pygna.command.read_distance_matrix
            import pygna.reading_class as rc
            import pygna.centrality.calculate_centrality

            matrix_id = {"nodes": read_distance_matrix(rwr_matrix_filename, in_memory=in_memory)[0],
               "matrix": read_distance_matrix(rwr_matrix_filename, in_memory=in_memory)[1]}
            graph = rc.ReadTsv(network_file).get_network()

            graph_centrality = calculate_centrality(graph, matrix_id)
    """

    graph_centrality = {}
    for n in graph.nodes:
        matrix_id = matrix["nodes"].index(n)
        sp = matrix["matrix"][matrix_id]
        tot_sp = sum(sp)
        if tot_sp > 0:
            # Remove 1 because we are considering one node of the graph
            graph_centrality[n] = (len(sp) - 1) / tot_sp

    return graph_centrality
