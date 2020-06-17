Command Line Tool
======================

Network Summaries and Visualisation
-----------------------------------

Before running any geneset-network analysis is a good practice to extract basic information on the network and the geneset or to visualise the network.
We provide a function to obtain the network summary ( or the summary of a single genesets) and another utility to write an annotated graphml file ( visualise it on Cytoscape ).

.. autofunction:: pygna.command.network_summary


.. code-block:: text

    pygna network-summary [-h] [-g GENESET_INPUT_FILE] [-s SETNAME] network-file text-output degree-figure-file c-components-figure-file

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

    usage: pygna build-distance-matrix [-h] [-g] network-file output-file

    positional arguments:
      network-file          network file
      output-file           distance matrix output file, use .hdf5

    optional arguments:
      -h, --help            show this help message and exit
      -g, --giant-component-only
                            compute the shortest paths only for nodes in the giant component (default: True)


Diffusion matrix
++++++++++++++++

To evaluate the diffusion matrix we have the below function that
implements a Randowm Walk with Restart algorithm. The $beta$ parameter
is set to $0.80$ as default, but can be given by the user.

.. autofunction:: pygna.command.build_rwr_diffusion

.. code-block:: text

    usage: pygna build-rwr-diffusion [-h] [-b BETA] [-o OUTPUT_FILE] network-file

    positional arguments:
      network-file          network file

    optional arguments:
      -h, --help            show this help message and exit
      -b BETA, --beta BETA  0.85
      -o OUTPUT_FILE, --output-file OUTPUT_FILE
                            distance matrix output file (use .hdf5) (default: -)

Show Results
------------

.. autofunction:: pygna.painter.paint_datasets_stats

.. code-block:: text

    usage: pygna paint-datasets-stats [-h] [-a ALTERNATIVE] table-filename output-file

    positional arguments:
      table-filename        pygna results table
      output-file           figure file, use pdf or png extension

    optional arguments:
      -h, --help            show this help message and exit
      -a ALTERNATIVE, --alternative ALTERNATIVE
                            'greater'

.. autofunction:: pygna.painter.paint_comparison_matrix

.. code-block:: text

    usage: pygna paint-comparison-matrix [-h] [-r] [-s] [-a] table-filename output-file

    positional arguments:
      table-filename        pygna comparison output
      output-file           output figure file, specify png or pdf file

    optional arguments:
      -h, --help            show this help message and exit
      -r, --rwr             use rwr is the table comes from a rwr analysis (default: False)
      -s, --single-geneset  use true if the comparison has been done for a single file (default: False)
      -a, --annotate        set true if uou want to print the pvalue inside the cell (default: False)


Utilities
---------


.. autofunction:: pygna.utils.geneset_from_table

.. code-block:: text

    usage: pygna geneset-from-table [-h] [--output-gmt OUTPUT_GMT] [--output-csv OUTPUT_CSV] [-n NAME_COLUMN] [-f FILTER_COLUMN] [-a ALTERNATIVE] [-t THRESHOLD]
                                [-d DESCRIPTOR]
                                input-file setname

    positional arguments:
      input-file            input csv file
      setname               name of the set

    optional arguments:
      -h, --help            show this help message and exit
      --output-gmt OUTPUT_GMT
                            output gmt name (default: -)
      --output-csv OUTPUT_CSV
                            output csv name (default: -)
      -n NAME_COLUMN, --name-column NAME_COLUMN
                            column with the names (default: 'Unnamed: 0')
      -f FILTER_COLUMN, --filter-column FILTER_COLUMN
                            column with the values to be filtered (default: 'padj')
      -a ALTERNATIVE, --alternative ALTERNATIVE
                            alternative to use for the filter, with less the filter is applied <threshold, otherwise >= threshold (default: 'less')
      -t THRESHOLD, --threshold THRESHOLD
                            threshold for the filter (default: 0.01)
      -d DESCRIPTOR, --descriptor DESCRIPTOR
                            descriptor for the gmt file (default: -)


.. autofunction:: pygna.utils.convert_gmt

.. code-block:: text

    usage: pygna convert-gmt [-h] [-e ENTREZ_COL] [-s SYMBOL_COL] gmt-file output-gmt-file conversion converter-map-filename

    positional arguments:
      gmt-file              gmt file to be converted
      output-gmt-file       output file
      conversion            e2s or s2e
      converter-map-filename
                            tsv table used to convert gene names

    optional arguments:
      -h, --help            show this help message and exit
      -e ENTREZ_COL, --entrez-col ENTREZ_COL
                            name of the entrez column (default: 'NCBI Gene ID')
      -s SYMBOL_COL, --symbol-col SYMBOL_COL
                            name of the symbol column (default: 'Approved symbol')


.. autofunction:: pygna.utils.generate_group_gmt

.. code-block:: text

    usage: pygna generate-group-gmt [-h] [-n NAME_COL] [-g GROUP_COL] [-d DESCRIPTOR] input-table output-gmt

    positional arguments:
      input-table           table to get the geneset from
      output-gmt            output gmt file

    optional arguments:
      -h, --help            show this help message and exit
      -n NAME_COL, --name-col NAME_COL
                            'Gene'
      -g GROUP_COL, --group-col GROUP_COL
                            'Cancer'
      -d DESCRIPTOR, --descriptor DESCRIPTOR
                            'cancer_genes'

.. autofunction:: pygna.utils.convert_csv

.. code-block:: text

    usage: pygna convert-csv [-h] [--converter-map-filename CONVERTER_MAP_FILENAME] [--output-file OUTPUT_FILE] [-e ENTREZ_COL] [-s SYMBOL_COL]
                         csv-file conversion original-name-col new-name-col geneset

    positional arguments:
      csv-file              csv file where to add a name column
      conversion            e2s or s2e
      original-name-col     column name to be converted
      new-name-col          name of the new column with the converted names
      geneset               the geneset to convert

    optional arguments:
      -h, --help            show this help message and exit
      --converter-map-filename CONVERTER_MAP_FILENAME
                            tsv table used to convert gene names (default: 'entrez_name.tsv')
      --output-file OUTPUT_FILE
                            if none, table is saved in the same input file (default: -)
      -e ENTREZ_COL, --entrez-col ENTREZ_COL
                            name of the entrez column (default: 'NCBI Gene ID')
      -s SYMBOL_COL, --symbol-col SYMBOL_COL
                            name of the symbol column (default: 'Approved symbol')

Block Model
___________

.. autofunction:: pygna.block_model.generate_sbm2_network

.. code-block:: text

    usage: pygna generate-sbm2-network [-h] [--prefix PREFIX] [--n-nodes N_NODES] [-t THETA0] [--percentage PERCENTAGE] [-d DENSITY] [--n-simulations N_SIMULATIONS]
                                   output-folder

    positional arguments:
      output-folder         folder where the simulations are saved

    optional arguments:
      -h, --help            show this help message and exit
      --prefix PREFIX       prefix for the simulations (default: 'sbm')
      --n-nodes N_NODES     nodes in the network (default: 1000)
      -t THETA0, --theta0 THETA0
                            probability of connection in the cluster (default: '0.9,0.7,0.5,0.2')
      --percentage PERCENTAGE
                            percentage of nodes in cluster 0, use ratio 0.1 = 10 percent (default: '0.1')
      -d DENSITY, --density DENSITY
                            multiplicative parameter used to define network density (default: '0.06,0.1,0.2')
      --n-simulations N_SIMULATIONS
                            number of simulated networks for each configuration (default: 3)


Painter
_______


.. autofunction:: pygna.painter.paint_datasets_stats

.. code-block:: text

    usage: pygna paint-datasets-stats [-h] [-a ALTERNATIVE] table-filename output-file

    positional arguments:
      table-filename        pygna results table
      output-file           figure file, use pdf or png extension

    optional arguments:
      -h, --help            show this help message and exit
      -a ALTERNATIVE, --alternative ALTERNATIVE
                            'greater'

.. autofunction:: pygna.painter.paint_comparison_matrix

.. code-block:: text

    usage: pygna paint-comparison-matrix [-h] [-r] [-s] [-a] table-filename output-file

    positional arguments:
      table-filename        pygna comparison output
      output-file           output figure file, specify png or pdf file

    optional arguments:
      -h, --help            show this help message and exit
      -r, --rwr             use rwr is the table comes from a rwr analysis (default: False)
      -s, --single-geneset  use true if the comparison has been done for a single file (default: False)
      -a, --annotate        set true if uou want to print the pvalue inside the cell (default: False)

.. autofunction:: pygna.painter.paint_volcano_plot

.. code-block:: text

    usage: pygna paint-volcano-plot [-h] [-r] [-i ID_COL] [--threshold-x THRESHOLD_X] [--threshold-y THRESHOLD_Y] [-a] table-filename output-file

    positional arguments:
      table-filename        pygna comparison output
      output-file           output figure file, specify png or pdf file

    optional arguments:
      -h, --help            show this help message and exit
      -r, --rwr             use rwr is the table comes from a rwr analysis (default: False)
      -i ID_COL, --id-col ID_COL
                            'setname_B'
      --threshold-x THRESHOLD_X
                            0
      --threshold-y THRESHOLD_Y
                            2
      -a, --annotate        False


.. autofunction:: pygna.painter.paint_summary_gnt

.. code-block:: text

    usage: pygna paint-summary-gnt [-h] [-s SETNAME] [-t THRESHOLD] [-c COLUMN_FILTER] [--larger] [--less-tests LESS_TESTS] output-figure [input_tables [input_tables ...]]

    positional arguments:
      output-figure         output figure filename
      input_tables          -

    optional arguments:
      -h, --help            show this help message and exit
      -s SETNAME, --setname SETNAME
                            name of the dataset (default: -)
      -t THRESHOLD, --threshold THRESHOLD
                            Value to threshold the colors (default: 0.05)
      -c COLUMN_FILTER, --column-filter COLUMN_FILTER
                            column where the threshold is applied (default: 'empirical_pvalue')
      --larger              if True the threshold is set as lower limit (default: False)
      --less-tests LESS_TESTS
                            comma separated string of the tests that are significant if lower than expected, otherwise pass empty string (default: 'topology_sp')
