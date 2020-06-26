Tutorial
========

General Usage
--------------

Data Input/Output
+++++++++++++++++

PyGNA relies on two fundamental objects, a network and a geneset. Regardless of the
test, the genes in the set are mapped to a given network (usually a large reference
network as the BioGRID one) and some statistic is evaluated.

Networks are read as tab separated text files with each couple of nodes that have an edge
between them: :math:`node_a  node_b`

Genesets use the gmt format, where each genesets:
`name \tab descriptor \tab gene1 \tab gene2`

Each gmt file can have a single geneset or multiple of them, PyGNA is able to both analyse all of them
or to restrict the analysis to a single geneset ( specifying the setname of interest ).

All the results are stored in csv tables reporting all results and parameters of the analysis.
Those can be easily visualised using the plotting utilities that produce either pdf or png figures.

Large Matrices
+++++++++++++++

Computing the shortest path and the interaction probabilities between each couple of nodes
in the network can be very expensive and time consuming. Luckily, this is a step
independent from the genesets and the analysis.

For this reason we have implemented a separate step for evaluating and saving the shortest path
and rwr matrices.

.. code-block:: bash

    $ pygna build-rwr-diffusion barabasi.interactome.tsv --output-file interactome_RWR.hdf5

.. code-block:: bash

    $ pygna build-distance-matrix barabasi.interactome.tsv interactome_SP.hdf5



Analysis
+++++++++++++++++

PyGNA provides both topology tests on single genesets or association tests on couples of them.
run pygna -h to show all possible options or check the whole :ref:`CLI` documentation.

Here a simplified structure of the available tests:

- **topology**:
    - module
    - internal degree
    - total degree
    - shortest path
    - random walk

- **association**:
    - shortest path
    - random walk

The call for each test is quite similar, we check now the topology-rwr test function
to review the fundamental parameters

.. code-block:: bash

    $ pygna test-topology-rwr -h

Complete Analysis of one geneset
--------------------------------

In case you have your own geneset you can completely characterise it using PyGNA as follows (names of min_working_example).

Necessary INPUT: <network> and <geneset> and possibly <geneset_pathway> to run the association analysis.
OUTPUT: <network_sp>.hdf5, <network_rwr>.hdf5, <table_results_test>_<test>.csv

Generate the network matrices:

.. code-block:: bash

    $ pygna build-distance-matrix <network> <network_sp>.hdf5
    $ pygna build-rwr-diffusion <network> --output-file <network_rwr>.hdf5

Topology tests:

.. code-block:: bash

    $ pygna test-topology-module <network> <geneset> <table_results_test>_topology_module.csv --number-of-permutations 100 --cores 4
    $ pygna test-topology-rwr <network> <geneset> <network_rwr>.hdf5 <table_results_test>_topology_rwr.csv --number-of-permutations 100 --cores 4
    $ pygna test-topology-internal-degree <network> <geneset> <table_results_test>_topology_internal_degree.csv --number-of-permutations 100 --cores 4
    $ pygna test-topology-sp <network> <geneset> <network_sp>.hdf5 <table_results_test>_topology_sp.csv --number-of-permutations 100 --cores 4
    $ pygna test-topology-total-degree <network> <geneset> <table_results_test>_topology_total_degree.csv --number-of-permutations 100 --cores 4

Association tests:

.. code-block:: bash

    $ pygna test-association-sp <network> <geneset> <network_sp>.hdf5 <table_results_test>_association_sp.csv -B <geneset_pathways> --keep --number-of-permutations 100 --cores 4
    $ pygna test-association-rwr <network> <geneset> <network_rwr>.hdf5 <table_results_test>_association_rwr.csv -B <geneset_pathways> --keep --number-of-permutations 100 --cores 4

Analysis of a geneset from a table (e.g. DeSeq2)
------------------------------------------------

We understand that in many cases the genes one wants to analyse are in a table-like format.
Hence, we provide a function to create a gmt geneset from a table, with the possibility of
applying a filter to the data. You can even just use it to return a gmt with all the genes
in a column by applying a dummy filter.

**NOTE**: In case you would like to apply more filters, just use the output_csv instead of
gmt, this way the first filters would just cut the original data returning the same table
format.

Here is how to obtain a gmt file of the significant genes obtained by DeSeq2.
we are here using *diff_exp* as name for the output geneset and we are filtering for padj<0.01.

.. code-block:: bash

    $ pygna geneset-from-table <deseq>.csv diff_exp <deseq>.gmt --descriptor deseq


Here is how to tweak the default behaviour to filter any csv table.

The filter is applied using the values in the filter_column (for example pvalues) and cutting using the
alternative and threshold parameters to specify what the filter should be. Bare in mind the filter
is supposed to be applied to **numerical values**. The output gmt will have the gene-names in the <name_column>

.. code-block:: bash

    $ pygna geneset-from-table <filename>.csv <setname> <filename>.gmt --name-colum <gene_names_column> --filter-column <filter-col> <'less'> --threshold <th> --descriptor <descriptor string>


Pipelines
---------

The package is integrable in Snakemake pipelines. We already provide
some sample analyses in [snakemake workflow](https://github.com/stracquadaniolab/workflow-pygna).
but we encourage to explore all the functionalities of the package and to raise issues
for bugs and alternative functionalities you might need.


Converting data using Pygna
+++++++++++++++++++++++++++

One of the most important feature in pygna is the possibility to convert a file from a format to another.
In particular here we propose some examples on how to use the converter classes.


Converting into GMT files
_________________________
From a general table containing all the genes to analyse, we can call the following command in order to get a GMT file that is correctly read from pygna:

.. code-block:: bash

    $ pygna geneset-from-table gene_analysis.csv brca --output-gmt brca_analysis.gmt -f significant -d significant -n genes.Entrezid -t 0.5 -a greater

It is possible to use pygna to merge different setnames in a single gmt file through the function `generate-group-gmt`.
You can override the default parameters, to match the reading exactly on your table.

.. code-block:: bash

    $ pygna generate-group-gmt gene_analysis.csv setnames_gmt.gmt group-col Cancer_Setnames

If you want just to add a column corresponding to the EntrezID or the gene's symbol, you can use the following command:

.. code-block:: bash

    $ pygna convert-csv mygenefile.csv e2s original-col-name EntrezID new-name-col Symbols geneset brca


Showing the results
--------------------

Pygna prepares ready-to-publish plots of all the analysis results.

Here is an example of a barplot of the GNT rwr analysis for multiple genesets:


.. code-block:: bash

    $ pygna paint-datasets-stats pygna_gnt_results.csv pygna_gnt.png


Which produced a plot similar to the following:


.. image:: _static/barplot.png

For a complete list of the plots refer to :ref:`visualisation`

Adding GNT or GNA test statistics
-----------------------------------

Pygna can be easily extended to perform different test statistics.
Check :ref:`customization` for a full tutorial on how to do that.

.. toctree::
    :maxdepth: 1

    vignettes

It is possible to define custom functions outside the pygna package and use them in a stand-alone test.
The code below shows how it is possible to implement a function such as the closeness centrality of a graph, using Pygna.

.. code-block:: python

    import argh
    import logging
    import networkx as nx
    import pygna.reading_class as rc
    import pygna.output as out
    import pygna.statistical_test as st
    import pygna.painter as paint
    import pygna.diagnostic as diagnostic
    import numpy as np


    def calculate_centrality(graph: nx.Graph, matrix: dict) -> np.ndarray:

        graph_centrality = list()
        for n in graph.nodes:
            matrix_id = matrix["nodes"].index(n)
            sp = matrix["matrix"][matrix_id]
            tot_sp = sum(sp)
            if tot_sp > 0:
                graph_centrality[n] = (len(sp) - 1) / tot_sp

        return np.asarray(graph_centrality)


    def test_topology_centrality(
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

        logging.info("Evaluating the test topology total degree, please wait")
        network = rc.ReadTsv(network_file).get_network()
        geneset = rc.ReadGmt(geneset_file).get_geneset(setname)
        setnames = [key for key in geneset.keys()]
        # Generate output
        output1 = out.Output(network_file, output_table, "topology_total_degree", geneset_file, setnames)
        logging.info("Results file = " + output1.output_table_results)
        # Create table
        output1.create_st_table_empirical()
        st_test = st.StatisticalTest(calculate_centrality, network)
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


Diagnostic
+++++++++++++++++++

For the GNT tests we also provide some diagnostic tools.

**#TODO: Add example of distribution plot**
