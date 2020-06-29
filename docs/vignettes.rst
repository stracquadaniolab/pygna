.. _customization:

Customizing the statistic elaboration
-----------------------------------------

In the following section, it is shown how to use the pygna classes in order to customize the statistic elaboration


Adding a GNT test
++++++++++++++++++++++++++

In addition to the provided functions, the user can define and use any function.
Follow the example below.


Step 0: define your function
____________________________
Define your custom function that implements your statistical test.
It must always return a `float`, representing the statistical value calculated in the function.

For example here it is reported the calculation of the geneset total degree

.. code-block:: python

    def geneset_total_degree_statistic(network, geneset):
        degree = nx.degree(network)
        geneset = list(geneset)
        total = np.array([degree[g] for g in geneset])
        return np.average(total)

Step 1: import the essential libraries
______________________________________

Import Pygna and its libraries

    >>> from pygna import *

Step 2: loading your data
_________________________

Load the network and the geneset you want to use, and initialise the output table

    >>> network = pygna.reading_class.ReadTsv("mynetwork.tsv").get_network()
    >>> geneset = rc.ReadGmt(geneset_file).get_geneset(setname)

    >>> setnames = [key for key in geneset.keys()]
    >>> output1 = out.Output(network_file, output_table, "topology_total_degree", geneset_file, setnames)


Step 4: setup the statistical comparison
________________________________________
Pass your custom function to Pygna statistical class like:

    >>> st_test = pygna.statistical_test.StatisticalTest(geneset_total_degree_statistic, network)

Step 5: apply the function to all the genes and update the output table
_________________________________________________________________________________________________

    >>> for setname, item in geneset.items():
    >>>    if len(item) > 20:
    >>>        item = set(item)
    >>>        observed, pvalue, null_d, n_mapped, n_geneset = st_test.empirical_pvalue(item,
    ...                                                                                 max_iter=number_of_permutations,
    ...                                                                                 alternative="greater",
    ...                                                                                 cores=cores)
    >>>        output1.update_st_table_empirical(setname, n_mapped, n_geneset, number_of_permutations, observed,
    ...                                   pvalue, np.mean(null_d), np.var(null_d))

Step 6: save the results
________________________


    >>>     output1.close_temporary_table()



Adding a GNA test
++++++++++++++++++++++++++


In addition to the provided functions, the user can define and use any function.
Follow the example below.


Step 0: define your function
____________________________

Before preparing the test define your custom function that implements your statistical test.
It must always return a `float`, representing the statistical value calculated in the function.

For example here it is reported the calculation of comparison shortest path between two genesets.
The return value is the average shortest path between the two genesets.

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


Step 1: import the essential libraries
______________________________________
Import Pygna and its libraries.

    >>> from pygna import *

Step 2: loading your data
_________________________
Load the network and the geneset you want to use, if needed load an hdf5 table

    >>> rw_dict = {"nodes": read_distance_matrix(rwr_matrix_filename, in_memory=in_memory)[0],
    ...            "matrix": read_distance_matrix(rwr_matrix_filename, in_memory=in_memory)[1]}
    >>> network = pygna.reading_class.ReadTsv("mynetwork.tsv").get_network()
    >>> geneset_a = rc.ReadGmt(file_geneset_a).get_geneset(setname_a)


Step 4: setup the statistical comparison
________________________________________
Pass your custom function to Pygna statistical class like:

    >>> st_comparison = pygna.statistical_comparison.StatisticalComparison(comparison_shortest_path, network, n_proc=cores, diz=rw_dict)

Step 5: apply the function to all the genes
___________________________________________
You can now perform the statistical comparison among all pairs in the geneset and update the output table

    >>> setnames = [key for key in geneset_a.keys()]
    >>> for pair in itertools.combinations(setnames, 2):
    >>>    if len(set(geneset_a[pair[0]])) > size_cut and len(set(geneset_a[pair[1]])) > size_cut:
    >>>        n_overlaps = len(set(geneset_a[pair[0]]).intersection(set(geneset_a[pair[1]])))
    >>>        observed, pvalue, null_d, a_mapped, b_mapped = st_comparison.comparison_empirical_pvalue(
    ...        set(geneset_a[pair[0]]), set(geneset_a[pair[1]]), max_iter=number_of_permutations, keep=keep)
    >>>        output1.update_comparison_table_empirical(pair[0], pair[1], len(set(geneset_a[pair[0]])), a_mapped,
    ...                                                          len(set(geneset_a[pair[1]])), b_mapped, n_overlaps,
    ...                                                      number_of_permutations, observed, pvalue, np.mean(null_d),
    ...                                                      np.var(null_d))

Step 6: save the results
________________________
Save the results using the pygna output function

    >>>     output1.close_temporary_table()
