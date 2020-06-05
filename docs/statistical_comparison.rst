Statistical Comparison
======================

.. autoclass:: pygna.statistical_comparison.StatisticalComparison
    :members:

.. autofunction:: pygna.statistical_comparison.comparison_shortest_path

.. autofunction:: pygna.statistical_comparison.calculate_sum

.. autofunction:: pygna.statistical_comparison.comparison_random_walk

Customizing the statistics
++++++++++++++++++++++++++

In addition to the provided functions, the user can define and use any function.
Follow the example below.

Step 1
______
Import Pygna and its libraries
    >>> from pygna import *

Step 2
______
Load the dictionary and the network you want to use
    >>> rw_dict = {"nodes": read_distance_matrix(rwr_matrix_filename, in_memory=in_memory)[0],
    ...            "matrix": read_distance_matrix(rwr_matrix_filename, in_memory=in_memory)[1]}
    >>> network = pygna.reading_class.ReadTsv("mynetwork.tsv").get_network()

Step 3
______
Define your custom function that implements your statistical test.
It must always return a `float` value.

For example here it is reported the calculation of comparison shortest path between two genesets.

.. code-block:: python

    def comparison_shortest_path(network, genesetA, genesetB, diz):
        n = np.array([diz["nodes"].index(i) for i in genesetA])
        m = np.array([diz["nodes"].index(i) for i in genesetB])
        cum_sum = calculate_sum(n, m, diz)
        cum_sum += calculate_sum(m, n, diz)
        d_AB = cum_sum / float(len(genesetA) + len(genesetB))
        d_A = st.geneset_localisation_statistic(network, genesetA, diz)
        d_B = st.geneset_localisation_statistic(network, genesetB, diz)
        return d_AB - (d_A + d_B) / 2

Step 4
______
Pass your custom function to Pygna statistical class like:
    >>> st_comparison = pygna.statistical_comparison.StatisticalComparison(comparison_shortest_path, network, n_proc=cores, diz=rw_dict)
