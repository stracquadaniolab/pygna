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
between them: node_a  node_b

Genesets use the gmt format, where each genesets has a name \tab descriptor \tab gene1 \tab gene2...

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

>>> pygna build-rwr-diffusion barabasi.interactome.tsv --output-file interactome_RWR.hdf5

>>> pygna build-distance-matrix barabasi.interactome.tsv interactome_SP.hdf5


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

>>> pygna test-topology-rwr -h

.. code-block::

    usage: pygna test-topology-rwr [-h] [--setname SETNAME] [--size-cut SIZE_CUT]
                                [--number-of-permutations NUMBER_OF_PERMUTATIONS]
                                [-c CORES] [-i]
                                [--results-figure RESULTS_FIGURE]
                                [-d DIAGNOSTIC_NULL_FOLDER]
                                network-file geneset-file rwr-matrix-filename
                                output-table

            Performs the analysis of random walk probabilities. 
            Given the RW matrix ( either normal random walk or RW with restart),
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
    --size-cut SIZE_CUT   removes all genesets with a mapped length < size_cut
                            (default: 20)
    --number-of-permutations NUMBER_OF_PERMUTATIONS
                            number of permutations for computing the empirical
                            pvalue (default: 500)
    -c CORES, --cores CORES
                            Number of cores for the multiprocessing (default: 1)
    -i, --in-memory       set if you want the large matrix to be read in memory
                            (default: False)
    --results-figure RESULTS_FIGURE
                            barplot of results, use pdf or png extension (default:
                            -)
    -d DIAGNOSTIC_NULL_FOLDER, --diagnostic-null-folder DIAGNOSTIC_NULL_FOLDER
                            plot null distribution, pass the folder where all the
                            figures are going to be saved (one for each dataset)
                            (default: -)

Complete Analysis of one geneset
--------------------------------

In case you have your own geneset you can completely characterise it using PyGNA as follows (names of min_working_example).

Necessary INPUT: <network> and <geneset> and possibly <geneset_pathway> to run the association analysis. 
OUTPUT: <network_sp>.hdf5, <network_rwr>.hdf5, <table_results_test>_<test>.csv 

Generate the network matrices:

>>> pygna build-distance-matrix <network> <network_sp>.hdf5
>>> pygna build-rwr-diffusion <network> --output-file <network_rwr>.hdf5

Topology tests:

>>> pygna test-topology-module <network> <geneset> <table_results_test>_topology_module.csv --number-of-permutations 100 --cores 4

>>> pygna test-topology-rwr <network> <geneset> <network_rwr>.hdf5 <table_results_test>_topology_rwr.csv --number-of-permutations 100 --cores 4

>>> pygna test-topology-internal-degree <network> <geneset> <table_results_test>_topology_internal_degree.csv --number-of-permutations 100 --cores 4


>>> pygna test-topology-sp <network> <geneset> <network_sp>.hdf5 <table_results_test>_topology_sp.csv --number-of-permutations 100 --cores 4


>>> pygna test-topology-total-degree <network> <geneset> <table_results_test>_topology_total_degree.csv --number-of-permutations 100 --cores 4

Association tests:

>>> pygna test-association-sp <network> <geneset> <network_sp>.hdf5 <table_results_test>_association_sp.csv -B <geneset_pathways> --keep --number-of-permutations 100 --cores 4

>>> pygna test-association-rwr <network> <geneset> <network_rwr>.hdf5 <table_results_test>_association_rwr.csv -B <geneset_pathways> --keep --number-of-permutations 100 --cores 4

Analysis of a geneset from a table (e.g. DeSeq2)
------------------------------------------------







Pipelines
---------

The package is integrable in Snakemake pipelines. We already provide
some sample analyses in [snakemake workflow](https://github.com/stracquadaniolab/workflow-pygna).
but we encourage to explore all the functionalities of the package and to raise issues
for bugs and alternative functionalities you might need.



