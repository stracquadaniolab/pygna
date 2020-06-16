Command Line Tool
======================

Network Summaries and Visualisation
-----------------------------------

Before running any geneset-network analysis is a good practice to extract basic information
on the network and the geneset or to visualise the network. We provide a function
to obtain the network summary ( or the summary of a single genesets) and another
utility to write an annotated graphml file ( visualise it on Cytoscape ).

.. autofunction:: pygna.command.network_summary

.. code-block:: text

    pygna network-summary [-h] [-g GENESET_INPUT_FILE] [-s SETNAME] network-file text-output degree-figure-file c-components-figure-file

    This function saves the principal info of a graph:
    - network properties
    - degree distribution
    - connected components diagnostic

    If a geneset/setname is passed to the function, the properties of
    the subgraph are evaluated


    positional arguments:
      network-file          network file
      text-output           output text file for the summary
      degree-figure-file    pdf or png file for the degree distribution
      c-components-figure-file
                            pdf or png file for the connected components distribution

    optional arguments:
      -h, --help            show this help message and exit
      -g GENESET_INPUT_FILE, --geneset-input-file GENESET_INPUT_FILE
                            geneset file (default: -)
      -s SETNAME, --setname SETNAME
                            specify a single geneset (default: -)


.. autofunction:: pygna.command.network_graphml

.. code-block:: txt

    usage: pygna network-graphml [-h] [-s SETNAME] [--giant-component-only] [-m] network-file geneset-file output-file

    This function generates a graphml file with nodes annotation.
    Given a geneset, with k setnames, each node has k False/True
    annotations for each set.

    Warning: without minimal, this function saves the full network.
    The minimal graph saves only the nodes in the geneset and those that
    connect them with a shortest path.


    positional arguments:
      network-file          network file
      geneset-file          geneset file
      output-file           graphml file for network for visualisation

    optional arguments:
      -h, --help            show this help message and exit
      -s SETNAME, --setname SETNAME
                            set name (default: -)
      --giant-component-only
                            saves only the giant component of the network (default: True)
      -m, --minimal         saves only the minimal graph (default: False)


Geneset Network Topology Analysis
---------------------------------

Here are all the GNT analyses that pygna can perform.
For each of them a single geneset topology is tested with the specified test statistics.


GNT Module analysis
+++++++++++++++++++

.. autofunction:: pygna.command.test_topology_module

.. code-block:: txt

    usage: pygna test-topology-module [-h] [--setname SETNAME] [--size-cut SIZE_CUT] [--number-of-permutations NUMBER_OF_PERMUTATIONS] [-c CORES] [--output-lcc OUTPUT_LCC]
                                  [-r RESULTS_FIGURE] [-d DIAGNOSTIC_NULL_FOLDER]
                                  network-file geneset-file output-table

    Performs geneset network topology module analysis.

    It computes a p-value for the largest connected component
    of the geneset being bigger than the one expected by chance
    for a geneset of the same size.


    positional arguments:
      network-file          network file
      geneset-file          GMT geneset file
      output-table          output results table, use .csv extension

    optional arguments:
      -h, --help            show this help message and exit
      --setname SETNAME     Geneset to analyse (default: -)
      --size-cut SIZE_CUT   removes all genesets with a mapped length < size_cut (default: 20)
      --number-of-permutations NUMBER_OF_PERMUTATIONS
                            number of permutations for computing the empirical pvalue (default: 500)
      -c CORES, --cores CORES
                            Number of cores for the multiprocessing (default: 1)
      --output-lcc OUTPUT_LCC
                            for creating a GMT file with the LCC lists pass a gmt filename (default: -)
      -r RESULTS_FIGURE, --results-figure RESULTS_FIGURE
                            barplot of results, use pdf or png extension (default: -)
      -d DIAGNOSTIC_NULL_FOLDER, --diagnostic-null-folder DIAGNOSTIC_NULL_FOLDER
                            plot null distribution, pass the folder where all the figures are going to be saved (one for each dataset) (default: -)


GNT Degree analysis
+++++++++++++++++++

.. autofunction:: pygna.command.test_topology_internal_degree

.. code-block:: txt

    usage: pygna test-topology-internal-degree [-h] [--setname SETNAME] [--size-cut SIZE_CUT] [--number-of-permutations NUMBER_OF_PERMUTATIONS] [-c CORES] [-r RESULTS_FIGURE]
                                           [-d DIAGNOSTIC_NULL_FOLDER]
                                           network-file geneset-file output-table

        Performs the analysis of internal degree.
        It computes a p-value for the ratio of internal degree
        of the geneset being bigger than the one expected by chance
        for a geneset of the same size.


    positional arguments:
      network-file          network file
      geneset-file          GMT geneset file
      output-table          output results table, use .csv extension

    optional arguments:
      -h, --help            show this help message and exit
      --setname SETNAME     Geneset to analyse (default: -)
      --size-cut SIZE_CUT   removes all genesets with a mapped length < size_cut (default: 20)
      --number-of-permutations NUMBER_OF_PERMUTATIONS
                            number of permutations for computing the empirical pvalue (default: 500)
      -c CORES, --cores CORES
                            Number of cores for the multiprocessing (default: 1)
      -r RESULTS_FIGURE, --results-figure RESULTS_FIGURE
                            barplot of results, use pdf or png extension (default: -)
      -d DIAGNOSTIC_NULL_FOLDER, --diagnostic-null-folder DIAGNOSTIC_NULL_FOLDER
                            plot null distribution, pass the folder where all the figures are going to be saved (one for each dataset) (default: -)

.. autofunction:: pygna.command.test_topology_total_degree

.. code-block:: txt

    usage: pygna test-topology-total-degree [-h] [--setname SETNAME] [--size-cut SIZE_CUT] [--number-of-permutations NUMBER_OF_PERMUTATIONS] [-c CORES] [-r RESULTS_FIGURE]
                                        [-d DIAGNOSTIC_NULL_FOLDER]
                                        network-file geneset-file output-table

        Performs the analysis of total degree of the .

        It computes a p-value for the ratio of total degree
        of the geneset being bigger than the one expected by chance
        for a geneset of the same size.


    positional arguments:
      network-file          network file
      geneset-file          GMT geneset file
      output-table          output results table, use .csv extension

    optional arguments:
      -h, --help            show this help message and exit
      --setname SETNAME     Geneset to analyse (default: -)
      --size-cut SIZE_CUT   removes all genesets with a mapped length < size_cut (default: 20)
      --number-of-permutations NUMBER_OF_PERMUTATIONS
                            number of permutations for computing the empirical pvalue (default: 500)
      -c CORES, --cores CORES
                            Number of cores for the multiprocessing (default: 1)
      -r RESULTS_FIGURE, --results-figure RESULTS_FIGURE
                            barplot of results, use pdf or png extension (default: -)
      -d DIAGNOSTIC_NULL_FOLDER, --diagnostic-null-folder DIAGNOSTIC_NULL_FOLDER
                        plot null distribution, pass the folder where all the figures are going to be saved (one for each dataset) (default: -)

GNT Shortest Path Analysis
+++++++++++++++++++++++++++

.. autofunction:: pygna.command.test_topology_sp

.. code-block:: text

    usage: pygna test-topology-sp [-h] [--setname SETNAME] [--size-cut SIZE_CUT] [--number-of-permutations NUMBER_OF_PERMUTATIONS] [-c CORES] [-i] [-r RESULTS_FIGURE]
                              [--diagnostic-null-folder DIAGNOSTIC_NULL_FOLDER]
                              network-file geneset-file distance-matrix-filename output-table

        Performs geneset network topology shortest path analysis.

        It computes a p-value for the average shortest path length
        of the geneset being smaller than expected by chance
        for a geneset of the same size.



    positional arguments:
      network-file          network file
      geneset-file          GMT geneset file
      distance-matrix-filename
                            distance hdf5 matrix file generated by pygna
      output-table          output results table, use .csv extension

    optional arguments:
      -h, --help            show this help message and exit
      --setname SETNAME     Geneset to analyse (default: -)
      --size-cut SIZE_CUT   removes all genesets with a mapped length < size_cut (default: 20)
      --number-of-permutations NUMBER_OF_PERMUTATIONS
                            number of permutations for computing the empirical pvalue (default: 500)
      -c CORES, --cores CORES
                            Number of cores for the multiprocessing (default: 1)
      -i, --in-memory       set if you want the large matrix to be read in memory (default: False)
      -r RESULTS_FIGURE, --results-figure RESULTS_FIGURE
                            barplot of results, use pdf or png extension (default: -)
      --diagnostic-null-folder DIAGNOSTIC_NULL_FOLDER
                            plot null distribution, pass the folder where all the figures are going to be saved (one for each dataset) (default: -)


GNT Random Walk Analysis
++++++++++++++++++++++++

.. autofunction:: pygna.command.test_topology_rwr

.. code-block:: text
    usage: pygna test-topology-rwr [-h] [--setname SETNAME] [--size-cut SIZE_CUT] [--number-of-permutations NUMBER_OF_PERMUTATIONS] [-c CORES] [-i]
                               [--results-figure RESULTS_FIGURE] [-d DIAGNOSTIC_NULL_FOLDER]
                               network-file geneset-file rwr-matrix-filename output-table

        Performs the analysis of random walk probabilities.
        Given the RWR matrix ,
        it compares the probability of walking between the genes in the geneset
        compared to those of walking between the nodes
        of a geneset with the same size


    positional arguments:
      network-file          network file, use a network with weights
      geneset-file          GMT geneset file
      rwr-matrix-filename   hdf5 RWR matrix obtained with pygna
      output-table          output results table, use .csv extension

    optional arguments:
      -h, --help            show this help message and exit
      --setname SETNAME     Geneset to analyse (default: -)
      --size-cut SIZE_CUT   removes all genesets with a mapped length < size_cut (default: 20)
      --number-of-permutations NUMBER_OF_PERMUTATIONS
                            number of permutations for computing the empirical pvalue (default: 500)
      -c CORES, --cores CORES
                            Number of cores for the multiprocessing (default: 1)
      -i, --in-memory       set if you want the large matrix to be read in memory (default: False)
      --results-figure RESULTS_FIGURE
                            barplot of results, use pdf or png extension (default: -)
      -d DIAGNOSTIC_NULL_FOLDER, --diagnostic-null-folder DIAGNOSTIC_NULL_FOLDER
                            plot null distribution, pass the folder where all the figures are going to be saved (one for each dataset) (default: -)


Geneset Network Association Analysis
------------------------------------

Here are all the GNA analyses that PyGNA can perform.
For each of them two sets are compared, one can pass:
- only one geneset: association between each set is computed
- two genesets: association between all the terms in the two genesets is computed

GNA Shortest Path
+++++++++++++++++++

.. autofunction:: pygna.command.test_association_sp

.. code-block:: text

    usage: pygna test-association-sp [-h] [--setname-a SETNAME_A] [--file-geneset-b FILE_GENESET_B] [--setname-b SETNAME_B] [--size-cut SIZE_CUT] [-k] [-c CORES] [-i]
                                 [--number-of-permutations NUMBER_OF_PERMUTATIONS] [-r RESULTS_FIGURE]
                                 network-file file-geneset-a distance-matrix-filename output-table

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


    positional arguments:
      network-file          network file
      file-geneset-a        GMT geneset file, if it's the only parameter passed the analysis is gonna be run on all the couples of datasets, otherwise specify the other files
                            and setnames
      distance-matrix-filename
                            distance matrix file generated by pygna
      output-table          output results table, use .csv extension

    optional arguments:
      -h, --help            show this help message and exit
      --setname-a SETNAME_A
                            Geneset A to analyse (default: -)
      --file-geneset-b FILE_GENESET_B
                            GMT geneset file (default: -)
      --setname-b SETNAME_B
                            Geneset B to analyse (default: -)
      --size-cut SIZE_CUT   removes all genesets with a mapped length < size_cut (default: 20)
      -k, --keep            if true, keeps the geneset B not permuted (default: False)
      -c CORES, --cores CORES
                            Number of cores for the multiprocessing (default: 1)
      -i, --in-memory       set if you want the large matrix to be read in memory (default: False)
      --number-of-permutations NUMBER_OF_PERMUTATIONS
                            number of permutations for computing the empirical pvalue (default: 500)
      -r RESULTS_FIGURE, --results-figure RESULTS_FIGURE
                            barplot of results, use pdf or png extension (default: -)


GNA Random Walk
+++++++++++++++++++

.. autofunction:: pygna.command.test_association_rwr

.. code-block:: text

    usage: pygna test-association-rwr [-h] [--setname-a SETNAME_A] [--file-geneset-b FILE_GENESET_B] [--setname-b SETNAME_B] [--size-cut SIZE_CUT] [-k] [-c CORES] [-i]
                                  [--number-of-permutations NUMBER_OF_PERMUTATIONS] [--results-figure RESULTS_FIGURE]
                                  network-file file-geneset-a rwr-matrix-filename output-table

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


    positional arguments:
      network-file          network file
      file-geneset-a        GMT geneset file
      rwr-matrix-filename   .hdf5 file with the RWR matrix obtained by pygna
      output-table          output results table, use .csv extension

    optional arguments:
      -h, --help            show this help message and exit
      --setname-a SETNAME_A
                            Geneset A to analyse (default: -)
      --file-geneset-b FILE_GENESET_B
                            GMT geneset file (default: -)
      --setname-b SETNAME_B
                            Geneset B to analyse (default: -)
      --size-cut SIZE_CUT   removes all genesets with a mapped length < size_cut (default: 20)
      -k, --keep            if true, keeps the geneset B unpermuted (default: False)
      -c CORES, --cores CORES
                            Number of cores for the multiprocessing (default: 1)
      -i, --in-memory       set if you want the large matrix to be read in memory (default: False)
      --number-of-permutations NUMBER_OF_PERMUTATIONS
                            number of permutations for computing the empirical pvalue (default: 500)
      --results-figure RESULTS_FIGURE
                            heatmap of results (default: -)

Weights Diffusion Analysis
--------------------------

**Please be aware that this function is still a beta version.**

.. autofunction:: pygna.command.test_diffusion_hotnet

.. code-block:: text

    usage: pygna test-diffusion-hotnet [-h] [--name-column NAME_COLUMN] [-w WEIGHT_COLUMN] [--filter-column FILTER_COLUMN] [--filter-condition FILTER_CONDITION]
                                   [--filter-threshold FILTER_THRESHOLD] [--normalise] [-s SIZE_CUT] [--number-of-permutations NUMBER_OF_PERMUTATIONS] [-c CORES] [-i]
                                   network-file geneset-file rwr-matrix-filename output-table

        Performs the analysis of random walk applying the weights of an upstream analysis.
        Given a csv file the user needs to specify the columns of interest and
        the threshold of significance.
        For the analysis the StatisticalDiffusion is used with hotnet_diffusion_statistic
        function.


    positional arguments:
      network-file          network file, use a network with weights
      geneset-file          csv geneset file
      rwr-matrix-filename   hdf5 RWR matrix obtained with pygna
      output-table          output results table, use .csv extension

    optional arguments:
      -h, --help            show this help message and exit
      --name-column NAME_COLUMN
                            Column to use as name (default is deseq2) (default: 'gene_name')
      -w WEIGHT_COLUMN, --weight-column WEIGHT_COLUMN
                            Column to use as weight (default is deseq2) (default: 'stat')
      --filter-column FILTER_COLUMN
                            Column used to define the significant genes (default is deseq2) (default: 'padj')
      --filter-condition FILTER_CONDITION
                            Condition for significance (default: 'less')
      --filter-threshold FILTER_THRESHOLD
                            threshold for significance (default: 0.01)
      --normalise           pass this flag for using only positive values in the analysis (default: False)
      -s SIZE_CUT, --size-cut SIZE_CUT
                            removes all genesets with a mapped length < size_cut (default: 20)
      --number-of-permutations NUMBER_OF_PERMUTATIONS
                            number of permutations for computing the empirical pvalue (default: 500)
      -c CORES, --cores CORES
                            Number of cores for the multiprocessing (default: 1)
      -i, --in-memory       set if you want the large matrix to be read in memory (default: False)


Large Matrices creation
-----------------------

Along with the network and the geneset, some analyses require an additional large matrix to be passed as input.
In particular, the analysis of shortest path and diffusion are evaluating a matrix of shape NxN (N being the number of nodes),
since those are invariant to the geneset analysed, they must be evaluated and saved only once.


Shortest Paths matrix
+++++++++++++++++++++

To evaluate the matrix of shortest paths we can use the following function:

.. autofunction:: pygna.command.build_distance_matrix

.. code-block:: text

Diffusion matrix
++++++++++++++++

To evaluate the diffusion matrix we have the below function that
implements a Randowm Walk with Restart algorithm. The $beta$ parameter
is set to $0.80$ as default, but can be given by the user.

.. autofunction:: pygna.command.build_rwr_diffusion

Show Results
------------

.. autofunction:: pygna.painter.paint_datasets_stats
.. autofunction:: pygna.painter.paint_comparison_matrix

Utilities
---------

.. autofunction:: pygna.utils.convert_gmt
.. autofunction:: pygna.utils.geneset_from_table
.. autofunction:: pygna.utils.filter_table
.. autofunction:: pygna.utils.generate_group_gmt
.. autofunction:: pygna.utils.convert_csv


Block Model
___________

.. autofunction:: pygna.block_model.generate_sbm2_network


Painter
_______


.. autofunction:: pygna.painter.paint_datasets_stats
.. autofunction:: pygna.painter.paint_comparison_matrix
.. autofunction:: pygna.painter.paint_volcano_plot
.. autofunction:: pygna.painter.plot_adjacency
