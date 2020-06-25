Command Line Tool
======================

Network Summaries and Visualisation
-----------------------------------

Before running any geneset-network analysis is a good practice to extract basic information on the network and the geneset or to visualise the network.
We provide a function to obtain the network summary ( or the summary of a single genesets) and another utility to write an annotated graphml file ( visualise it on Cytoscape ).

Network Summary
***************

    This function saves the principal info of a graph: - network properties - degree distribution - connected components diagnostic

    If a geneset/setname is passed to the function, the properties of the subgraph are evaluated

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


Network Graphml
***************

    This function generates a graphml file with nodes annotation. Given a geneset, with k setnames, each node has k False/True annotations for each set.

    Warning: without minimal, this function saves the full network. The minimal graph saves only the nodes in the geneset and those that connect them with a shortest path.

.. code-block:: text

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

Test Topology Module
********************

    Performs geneset network topology module analysis.

    It computes a p-value for the largest connected component of the geneset being bigger than the one expected by chance for a geneset of the same size.


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

Test topology internal degree
*****************************


Performs the analysis of internal degree. It computes a p-value for the ratio of internal degree of the geneset being bigger than the one expected by chance for a geneset of the same size.


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

Test topology total degree
**************************

    Performs the analysis of total degree of the network.

    It computes a p-value for the ratio of total degree of the geneset being bigger than the one expected by chance for a geneset of the same size.


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

Test topology shortest path
***************************

    Performs geneset network topology shortest path analysis.

    It computes a p-value for the average shortest path length of the geneset being smaller than expected by chance for a geneset of the same size.


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

Test topology random walk with restart
**************************************


Performs the analysis of random walk probabilities. Given the RWR matrix, it compares the probability of walking between the genes in the geneset compared to those of walking between the nodes of a geneset with the same size


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

Test association shortest path
******************************

    Performs comparison of network location analysis. If the flag –keep is passed, the B geneset is kept fixed, and doesnt’t get permuted

    It computes a p-value for the shortest path distance between two genesets being smaller than expected by chance

    If only A_geneset_file is passed the analysis is run on all the couples of sets in the file, if both A_geneset_file and B_geneset_file are passed, one can specify the setnames for both, if there is only one geneset in the file, setname_X can be omitted, if both sets are in the same file, B_geneset_file can be not specified, but setnames are needed.


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

Test association random walk with restart
*****************************************


    Performs comparison of network location analysis.

    It computes a p-value for the shortest path distance between two genesets being smaller than expected by chance

    If only A_geneset_file is passed the analysis is run on all the couples of sets in the file, if both A_geneset_file and B_geneset_file are passed, one can specify the setnames for both, if there is only one geneset in the file, setname_X can be omitted, if both sets are in the same file, B_geneset_file can be not specified, but setnames are needed.


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

Test diffusion hotnet
*********************

Performs the analysis of random walk applying the weights of an upstream analysis. Given a csv file the user needs to specify the columns of interest and the threshold of significance. For the analysis the StatisticalDiffusion is used with hotnet_diffusion_statistic function.


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

Build distance matrix
*********************

Build a shortest path distance matrix for a given network. Matrix can be saved as a .txt file or a .hdf5 one.

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

Build random walk diffusion matrix
**********************************

Build the RWR_diffusion_matrix

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

Paint dataset statistics
************************

This function plots the results of of a GNT test. Pass the results table generated by one of the functions and the output figure file (png or pdf). In case you are using a SP test, pass also ‘less’ as an alternative.


.. code-block:: text

    usage: pygna paint-datasets-stats [-h] [-a ALTERNATIVE] table-filename output-file

    positional arguments:
      table-filename        pygna results table
      output-file           figure file, use pdf or png extension

    optional arguments:
      -h, --help            show this help message and exit
      -a ALTERNATIVE, --alternative ALTERNATIVE
                            'greater'

Paint comparison matrix
***********************

This function plots the results of of a GNA test. Pass the results table generated by one of the functions and the output figure file (png or pdf). With rwr you can specify whether the test is a rwr association, in this case a different palette and limits are sets. Specify if the results are obtained using association with only one genesets (multiple setnames in the same file). Pass the annotate flag to have the pvalue annotation on the plot

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

Geneset from table
******************

Generate a table with the genes setname.

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


Convert gmt file
****************

Convert a GMT file by adding a the entrez or symbol column for each gene

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


Generate group of gmt
*********************

Generate a single GMT file with the addition of the column description

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

Convert csv file
****************

Convert a CSV file by adding a conversion genes column (for entrezID or symbol)

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

Get connected components
************************

This function evaluate all the connected components in the subgraph pf the network with a given setname. Multiple setnames can be passed to this function to analyze all of them in a run. The file produces a gmt output and optionally a plot of the subnetwork with the connected components analysed.

Please notice that to convert the entrezID into Symbols, a stable internet connection is required


.. code-block:: text

    usage: pygna get-connected-components [-h] [-c] network-file geneset-file setname o graphml

    positional arguments:
      network-file          network file
      geneset-file          GMT geneset file
      setname               The setname to analyse
      o                     The output file name (should be gmt)
      graphml               The name of the graphml file

    optional arguments:
      -h, --help            show this help message and exit
      -c, --convert-entrez  Convert EntrezID->Symbol (default: True)


Block Model
___________

Generate stochastic block model network
***************************************

This function generates the simulated networks and genesets using the stochastic block model with 2 BLOCKS as described in the paper. The output names are going to be prefix_t_<theta0>_p_<percentage>_d_<density>_s_<n_simulation>_network.tsv or _genes.gmt One connected cluster while the rest of the network has the same probability of connection. SBM = d *[theta0, 1-theta0 1-theta0, 1-theta0] The simulator checks for connectedness of the generated network, if the generated net is not connected, a new simulation is generated.

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

Paint comparison matrix
***********************

This function plots the results of of a GNA test. Pass the results table generated by one of the functions and the output figure file (png or pdf). With rwr you can specify whether the test is a rwr association, in this case a different palette and limits are sets. Specify if the results are obtained using association with only one genesets (multiple setnames in the same file). Pass the annotate flag to have the pvalue annotation on the plot

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

Paint a volcano plot
********************

This function plots the results of of a GNA test of association of a single geneset against multiple pathways. Pass the results table generated by one of the functions and the output figure file (png or pdf). From the results table, a multiple testing correction is applied and the results are those plotted.

The defined threshold are for x: zscore and y: -log10(pvalue)


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


Paint summary of GNT analysis
*****************************

Plot the summary for the GNT analysis

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
