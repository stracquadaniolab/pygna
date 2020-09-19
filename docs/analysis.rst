Geneset Network Topology Analysis
---------------------------------

Here are all the GNT analyses that pygna can perform. For each of them a single
geneset topology is tested with the specified test statistics.


GNT Module analysis
+++++++++++++++++++

Test Topology Module
********************

`test-topology-module` performs geneset network topology module analysis.

It computes a p-value for the largest connected component of the geneset being bigger than the one expected by chance for a geneset of the same size.


.. code-block:: text

    usage: pygna test-topology-module [-h] [--setname SETNAME] [--size-cut SIZE_CUT] [--number-of-permutations NUMBER_OF_PERMUTATIONS] [-c CORES] [--n-bins N_BINS] [--output-lcc OUTPUT_LCC]
                                    [-r RESULTS_FIGURE] [-d DIAGNOSTIC_NULL_FOLDER]
                                    network-file geneset-file output-table

        Performs geneset network topology module analysis.
        It computes a p-value for the largest connected component of the geneset being bigger than the one expected by chance
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
    --n-bins N_BINS       if >1 applies degree correction by binning the node degrees and sampling according to geneset distribution (default: 1)
    --output-lcc OUTPUT_LCC
                            for creating a GMT file with the LCC lists pass a GMT filename (default: -)
    -r RESULTS_FIGURE, --results-figure RESULTS_FIGURE
                            barplot of results, use pdf or png extension (default: -)
    -d DIAGNOSTIC_NULL_FOLDER, --diagnostic-null-folder DIAGNOSTIC_NULL_FOLDER
                            plot null distribution, pass the folder where all the figures are going to be saved (one for each dataset) (default: -)


GNT Degree analysis
+++++++++++++++++++

Test topology internal degree
*****************************


`test-topology-internal-degree` performs the analysis of internal degree.
It computes a p-value for the ratio of internal degree of the geneset being
bigger than the one expected by chance for a geneset of the same size.

.. code-block:: text

    usage: pygna test-topology-internal-degree [-h] [--setname SETNAME] [--size-cut SIZE_CUT] [--number-of-permutations NUMBER_OF_PERMUTATIONS] [-c CORES] [--n-bins N_BINS] [-r RESULTS_FIGURE]
                                            [-d DIAGNOSTIC_NULL_FOLDER]
                                            network-file geneset-file output-table

        Performs the analysis of internal degree.
        It computes a p-value for the ratio of internal degree of the geneset being bigger than the one expected by chance
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
    --n-bins N_BINS       if >1 applies degree correction by binning the node degrees and sampling according to geneset distribution (default: 1)
    -r RESULTS_FIGURE, --results-figure RESULTS_FIGURE
                            barplot of results, use pdf or png extension (default: -)
    -d DIAGNOSTIC_NULL_FOLDER, --diagnostic-null-folder DIAGNOSTIC_NULL_FOLDER
                            plot null distribution, pass the folder where all the figures are going to be saved (one for each dataset) (default: -)


Test topology total degree
**************************

`test-topology-total-degree` performs the analysis of total degree of the network.
It computes a p-value for the ratio of total degree of the geneset being bigger than
the one expected by chance for a geneset of the same size.


.. code-block:: text

    usage: pygna test-topology-total-degree [-h] [--setname SETNAME] [--size-cut SIZE_CUT] [--number-of-permutations NUMBER_OF_PERMUTATIONS] [-c CORES] [-r RESULTS_FIGURE]
                                            [-d DIAGNOSTIC_NULL_FOLDER]
                                            network-file geneset-file output-table

        Performs the analysis of total degree of the .

        It computes a p-value for the ratio of total degree of the geneset being bigger than the one expected by chance
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
                            plot null distribution, pass the folder where all the figures
                            are going to be saved (one for each dataset) (default: -)


GNT Shortest Path Analysis
+++++++++++++++++++++++++++

Test topology shortest path
***************************

`test-topology-sp` performs geneset network topology shortest path analysis.
It computes a p-value for the average shortest path length of the geneset being smaller
than expected by chance for a geneset of the same size.

Consider that the computation of the shortest path matrix for a network can be considerably slow,
however it only needs to be done once for each network (and can then be reused).
We have separated the generation of the matrix step, refer to `build-distance-matrix`.


.. code-block:: text

    usage: pygna test-topology-sp [-h] [--setname SETNAME] [--size-cut SIZE_CUT] [--number-of-permutations NUMBER_OF_PERMUTATIONS] [-c CORES] [-i] [--n-bins N_BINS] [-r RESULTS_FIGURE]
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
    --n-bins N_BINS       if >1 applies degree correction by binning the node degrees and sampling according to geneset distribution (default: 1)
    -r RESULTS_FIGURE, --results-figure RESULTS_FIGURE
                            barplot of results, use pdf or png extension (default: -)
    --diagnostic-null-folder DIAGNOSTIC_NULL_FOLDER
                            plot null distribution, pass the folder where all the figures are going to be saved (one for each dataset) (default: -)


GNT Random Walk Analysis
++++++++++++++++++++++++

Test topology random walk with restart
**************************************


`test-topology-rwr` performs the analysis of random walk probabilities.
Given the Random Walk with Restart matrix, it compares the probability of walking between the genes
in the geneset compared to those of walking between the nodes of a geneset with the same size

For the generation of the RWR matrix check `build-rwr-diffusion`

.. code-block:: text

    usage: pygna test-topology-rwr [-h] [--setname SETNAME] [--size-cut SIZE_CUT] [--number-of-permutations NUMBER_OF_PERMUTATIONS] [-c CORES] [-i] [--n-bins N_BINS]
                                [--results-figure RESULTS_FIGURE] [-d DIAGNOSTIC_NULL_FOLDER]
                                network-file geneset-file rwr-matrix-filename output-table

        Performs the analysis of random walk probabilities.
        Given the RWR matrix, it compares the probability of walking between the genes in the geneset compared to
        those of walking between the nodes of a geneset with the same size


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
    --n-bins N_BINS       if >1 applies degree correction by binning the node degrees and sampling according to geneset distribution (default: 1)
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

`pygna test-association-sp` performs comparison of network location analysis. If the flag –keep is passed, the B geneset is kept fixed, and doesnt’t get permuted
It computes a p-value for the shortest path distance between two genesets being smaller than expected by chance
If only A_geneset_file is passed the analysis is run on all the pair of sets in the file, if both A_geneset_file and B_geneset_file are passed, one can specify the setnames for both, if there is only one geneset in the file, setname_X can be omitted, if both sets are in the same file, B_geneset_file can be not specified, but setnames are needed.


.. code-block:: text

    usage: pygna test-association-sp [-h] [--setname-a SETNAME_A] [--file-geneset-b FILE_GENESET_B] [--setname-b SETNAME_B] [--size-cut SIZE_CUT] [-k] [-c CORES] [-i]
                                    [--number-of-permutations NUMBER_OF_PERMUTATIONS] [--n-bins N_BINS] [-r RESULTS_FIGURE]
                                    network-file file-geneset-a distance-matrix-filename output-table

        Performs comparison of network location analysis. If the flag
        --keep  is passed, the B geneset is kept fixed, and doesnt't get permuted.

        It computes a p-value for the shortest path distance between two genesets being smaller than expected by chance
        If only A_geneset_file is passed the analysis is run on all the pair of sets in the file, if both
        A_geneset_file and B_geneset_file are passed, one can specify the setnames for both, if there is only one
        geneset in the file, setname_X can be omitted, if both sets are in the same file, B_geneset_file can be not
        specified, but setnames are needed.


    positional arguments:
    network-file          network file
    file-geneset-a        GMT geneset file, if it's the only parameter passed the analysis is gonna be run on all the pair of datasets, otherwise specify the other files and setnames
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
    --n-bins N_BINS       if >1 applies degree correction by binning the node degrees and sampling according to geneset distribution (default: 1)
    -r RESULTS_FIGURE, --results-figure RESULTS_FIGURE
                            barplot of results, use pdf or png extension (default: -)


GNA Random Walk
+++++++++++++++++++

Test association random walk with restart
*****************************************


`test-association-rwr` performs comparison of network location analysis.
It computes a p-value for the shortest path distance between two genesets being smaller than expected by chance
If only A_geneset_file is passed the analysis is run on all the couples of sets in the file, if both A_geneset_file and B_geneset_file are passed, one can specify the setnames for both, if there is only one geneset in the file, setname_X can be omitted, if both sets are in the same file, B_geneset_file can be not specified, but setnames are needed.


.. code-block:: text

    usage: pygna test-association-rwr [-h] [--setname-a SETNAME_A] [--file-geneset-b FILE_GENESET_B] [--setname-b SETNAME_B] [--size-cut SIZE_CUT] [-k] [-c CORES] [-i]
                                    [--number-of-permutations NUMBER_OF_PERMUTATIONS] [--n-bins N_BINS] [--results-figure RESULTS_FIGURE]
                                    network-file file-geneset-a rwr-matrix-filename output-table

        Performs comparison of network location analysis.

        It computes a p-value for the shortest path distance
        between two genesets being smaller than expected by chance.

        If only A_geneset_file is passed the analysis is run on all the pair of sets in the file, if both
        A_geneset_file and B_geneset_file are passed, one can specify the setnames for both, if there is only one
        geneset in the file, setname_X can be omitted, if both sets are in the same file, B_geneset_file can be not
        specified, but setnames are needed.


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
    --n-bins N_BINS       if >1 applies degree correction by binning the node degrees and sampling according to geneset distribution (default: 1)
    --results-figure RESULTS_FIGURE
                            heatmap of results (default: -)


Weights Diffusion Analysis
-----------------------------

Test diffusion hotnet
++++++++++++++++++++++

Performs the analysis of random walk applying the weights of an upstream analysis.
Given a csv file the user needs to specify the columns of interest and the threshold of significance.
For the analysis the StatisticalDiffusion is used with hotnet_diffusion_statistic function.


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




Hypothesis testing
--------------------

For both GNT and GNA we carry out hypothesis testing where the null hypothesis
is that the statistic for the tested geneset is not more extreme than expected by change.
Since closed form definitions for all null distributions are difficult to
define, we compute the null distributions by randomly resampling the nodes
of the network conditioned on the geneset size.

For example, for a GNT on a geneset of size :math:`n`
by passing the Number Of Permutations :math:`NOP` parameter,
:math:`NOP` random sets or nodes are drawn and on each of the test
statistic to build the null distribution.

Let :math:`\bar{q}` the observed value for the statistic and
:math:`Q` the null distribution we derive an empirical p-value as follows:

.. math::

    P(\bar{q} \geq Q )  = \dfrac{(\sum_{i=1}^{NOP}I(Q_i \geq \bar{q}))+1}{NOP+1}

where :math:`I` is the indicator function returning 1 if and only if the evaluated
condition is true and continuity correction is applied.


Degree corrected sampling
++++++++++++++++++++++++++++

We also allow to correct the resampling strategy on the degree distribution
of the input geneset. This way the null distribution is normalised
on the degree of the test set.

To enable this function pass the parameter `--n-bins Nb` with :math:`Nb > 1`.
In this case, the degrees of all nodes in the network are estimated and
their frequency computed on a binned histogram. The Nb parameter is indeed
used to specify the number of bins for the histogram. Consider that a low number
of bins corresponds to a low resolution of the histogram, however specifying too many
bins might correspond to sampling at random from the same nodes of the tested
geneset.

In particular the probability of sampling node :math:`y`
is proportional to (tested nodes in bin k)/((tested nodes) (all nodes in bin k))

GNA: association sampling
++++++++++++++++++++++++++++++

In a GNA two genesets are tested for their association.
When testing a signle geneset against many pathways we recommend
the use of the `--keep` flag. This way, while resampling only the geneset
a will be randomly permuted and the geneset b is going to be kept
as it is. This strategy is more conservative and is helpful
in testing whether the tested geneset is more strongly connected to
the pathway (or any other geneset of interest) than expected by chance.




