Tutorial
========

General Usage
--------------

Data Input/Output
+++++++++++++++++

PyGNA relies requires two data sources, a network and a geneset. Regardless of
the test, the genes in the geneset are mapped to the input network (e.g. BioGrid) and
then statistical tests are carried out.

Networks are read as tab separated text files (.tsv), where edges are represented by a node pair.
An example of a .tsv file is the fowllowing:

.. code-block:: text

    node1\tnode2
    node1\tnode3
    node1\tnode4
    node2\tnode5


Genesets use the GMT format, where each geneset is reported as:

.. code-block:: text

    <name>\t<descriptor>\tgene1\tgene2\tgene3...


Otherwise, a single geneset could be passed as a TXT file, with a list of genes, one for each line.

.. code-block:: text

    gene1
    gene2
    gene3
    ...

Since GMT file can have multiple genesets, PyGNA can either run the analyses on all
of them or on a user-specified subset.

Results are stored in CSV files, along with the parameters of the analysis.
Results can be easily visualised using PyGNA plotting utilities, and save either
as PDF or PNG files.

Matrix computation
++++++++++++++++++

Computing the shortest path and the interaction probabilities between each pair
of nodes in the network can be computationally taxing. However, since matrices
are never modified by statistical tests, they can be computed as part of a
pre-processing step.

For this reason, we have implemented a separate step for evaluating and saving
the shortest path and RWR matrices.

.. code-block:: bash

    $ pygna build-rwr-diffusion barabasi.interactome.tsv --output-file interactome_RWR.hdf5

.. code-block:: bash

    $ pygna build-distance-matrix barabasi.interactome.tsv interactome_SP.hdf5



Analysis
+++++++++++++++++

PyGNA provides commands for running geneset network topology (GNT) and geneset
network analysis (GNA) tests. Running `pygna -h` shows all the available
analyses, whereas extend documentation can be found :ref:`CLI` here.

Here a simplified structure of the available tests:

- **geneset network topology**:
    - module
    - internal degree
    - total degree
    - shortest path
    - random walk

- **geneset network association**:
    - shortest path
    - random walk

The analyses commands have all the same interface; for example, the available RWR GNT analysis
options can be checked by running:

.. code-block:: bash

    $ pygna test-topology-rwr -h

Complete analysis of one geneset
--------------------------------

In case you have your own geneset you can completely characterise it using PyGNA
as follows (names of min_working_example).

INPUT: <network> and <geneset>, and/or a <geneset_pathway> to run the
association analysis.

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

In many workflows, the genes to analyse are stored in a table-like format. Hence,
we provide a function to create a GMT geneset from a table, with the possibility
of applying a filter to the data. You can even just use it to return a GMT with
all the genes in a column by applying a dummy filter.

**NOTE**: In case you would like to apply more filters, just use the output_csv
command, instead of GMT, in order to only filter the original data and return
the table in the original format.

Here, for example, we obtain a GMT file of the differentially expressed genes
computed using DeSeq2, by selecting genes with padj<0.01.

.. code-block:: bash

    $ pygna geneset-from-table <deseq>.csv diff_exp <deseq>.gmt --descriptor deseq


PyGNA implements a generic interface to filter any CSV file. Filters are applied
to the values in the filter_column (for example pvalue) and the cut is performed using the
alternative parameter. The threshold parameters are used to specify what the filter should be used. Bare
in mind the filter is supposed to be applied to **numerical values**. The output
GMT will have the gene-names in the <name_column>

.. code-block:: bash

    $ pygna geneset-from-table <filename>.csv <setname> <filename>.gmt --name-colum <gene_names_column> --filter-column <filter-col> <'less'> --threshold <th> --descriptor <descriptor string>


Pipelines
---------

PyGNA can be seamlessly integrated into Snakemake workflows, and we provide a
basic [snakemake workflow](https://github.com/stracquadaniolab/workflow-pygna)
to carry out network analysis with PyGNA.

Converting data using PyGNA
+++++++++++++++++++++++++++

One of the most important feature in pygna is the possibility to convert a file
from a format to another. PyGNA supports:

Converting into GMT format
__________________________

Geneset in table-like formats (e.g. CSV) can be converted into GMT format as
follows:

.. code-block:: bash

    $ pygna geneset-from-table gene_analysis.csv brca --output-gmt brca_analysis.gmt -f significant -d significant -n genes.Entrezid -t 0.5 -a greater

It is also possible to merge different setnames in a single gmt file through the
function `generate-group-gmt`. You can override the default parameters, to match
the columns in your table.

.. code-block:: bash

    $ pygna generate-group-gmt gene_analysis.csv setnames_gmt.gmt group-col Cancer_Setnames

If you want to add a column corresponding to the EntrezID or the gene's symbol,
you can use the following command:

.. code-block:: bash

    $ pygna convert-csv mygenefile.csv e2s original-col-name EntrezID new-name-col Symbols geneset brca


Showing the results
--------------------

PyGNA generates publication-ready plots for each analysis.

For example, a barplot of the GNT RWR analysis for multiple genesets can be
generated by running:


.. code-block:: bash

    $ pygna paint-datasets-stats pygna_gnt_results.csv pygna_gnt.png


which produced a plot similar to the following:

.. image:: _static/barplot.png

For a complete list of the plots refer to :ref:`visualisation`

Adding GNT or GNA test statistics
-----------------------------------

PyGNA can be easily extended to run user-defined statistical tests. Check
" :ref:`customization` " for a full tutorial on how to do that.

.. toctree::
    :maxdepth: 1

    vignettes

It is possible to define custom functions outside the pygna package and use them in a stand-alone test.
The code below shows how it is possible to implement a function such as the average closeness centrality of a geneset, using Pygna.

.. code-block:: python

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
            This function calculates the average closeness centrality of a geneset.
            For a single node, the closeness centrality is defined as the inverse
            of the shortest path distance of the node from all the other nodes.

            Given a network with N nodes and a distance shortest path function
            between two nodes d(u,v)
            closeness centrality (u)= (N -1) / sum (v != u) d(u,v)

            where sp is the distance of the node with each other node and tot_sp is the total shortest paths for the whole graph.

            """

            graph_centrality = []

            ids = [diz["nodes"].index(n) for n in geneset]
            graph_centrality = [(len(matrix["nodes"]) - 1) / matrix['vector'][idx] for idx in ids]

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
            results_figure: "barplot of results, use pdf or png extension" = None,
            in_memory: 'load hdf5 data onto memory' = False,
            diagnostic_null_folder: "plot null distribution, pass the folder where all the figures are going to be saved "
                                    "(one for each dataset)" = None):

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
            output1 = out.Output(network_file, output_table, "topology_total_degree", geneset_file, setnames)
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





Diagnostic
---------------

Distribution plot
++++++++++++++++++++

When running a statistical test, one might want to visually assess the null distribution.
By passing  `-d <diagnostic_folder/>` through command line, a distribution plot
of the empirical null is shown for each test.

Here is an example.

.. image:: _static/distribution.png
