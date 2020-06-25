Geneset Network Topology Analysis
---------------------------------

Here are all the GNT analyses that pygna can perform.
For each of them a single geneset topology is tested with the specified test statistics.


GNT Module analysis
+++++++++++++++++++

.. autofunction:: pygna.command.test_topology_module

.. code-block:: text

    usage: pygna test-topology-module [-h] [--setname SETNAME] [--size-cut SIZE_CUT] [--number-of-permutations NUMBER_OF_PERMUTATIONS] [-c CORES] [--output-lcc OUTPUT_LCC]
                                  [-r RESULTS_FIGURE] [-d DIAGNOSTIC_NULL_FOLDER]
                                  network-file geneset-file output-table

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

.. code-block:: text

    usage: pygna test-topology-internal-degree [-h] [--setname SETNAME] [--size-cut SIZE_CUT] [--number-of-permutations NUMBER_OF_PERMUTATIONS] [-c CORES] [-r RESULTS_FIGURE]
                                           [-d DIAGNOSTIC_NULL_FOLDER]
                                           network-file geneset-file output-table

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

.. code-block:: text

    usage: pygna test-topology-total-degree [-h] [--setname SETNAME] [--size-cut SIZE_CUT] [--number-of-permutations NUMBER_OF_PERMUTATIONS] [-c CORES] [-r RESULTS_FIGURE]
                                        [-d DIAGNOSTIC_NULL_FOLDER]
                                        network-file geneset-file output-table

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


.. autofunction:: pygna.command.test_diffusion_hotnet

.. code-block:: text

    usage: pygna test-diffusion-hotnet [-h] [--name-column NAME_COLUMN] [-w WEIGHT_COLUMN] [--filter-column FILTER_COLUMN] [--filter-condition FILTER_CONDITION]
                                   [--filter-threshold FILTER_THRESHOLD] [--normalise] [-s SIZE_CUT] [--number-of-permutations NUMBER_OF_PERMUTATIONS] [-c CORES] [-i]
                                   network-file geneset-file rwr-matrix-filename output-table

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

