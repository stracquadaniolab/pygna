Quickstart
============

Installation
------------

The easiest and fastest way to install `pygna` using `conda`:

    >>> conda install -c stracquadaniolab -c bioconda -c conda-forge pygna

Alternatively you can install it through `pip`:

    >>> pip install pygna

Please note, that `pip` will not install non Python requirements.

Getting started
---------------

A typical `pygna` analysis consists of 3 steps:

1. Generate the RWR and SP matrices for the network you are using ( once they are generated, you won't need to repeat the same step again)
2. Make sure that the input genesets are in the right format. If a network uses entrez ID, and your file is in HUGO symbols, use the pygna utility for the name conversion.
3. Run the analysis you are interested into.
4. Once you have the output tables, you can choose to visualize one or more plots.

Otherwise you can check our [snakemake workflow](https://github.com/stracquadaniolab/workflow-pygna) for the full geneset analysis;
our workflow contains sample data that you can use to familiarize with our software.


The examples below show some basic analysis that can be carried out with pygna

Example 1: Running pygna GNT analysis
+++++++++++++++++++++++++++++++++++++

In the `min-working-example` we provide some basic files to run this minimal example. You can follow this instructions to run a network 
topology analysis on a single geneset file.
WARNING: *Pay attention to the fact we set the number of permutations to 1000, if you are just willing to test the behaviour, use 50 or 100
to speed up the process*

>>> cd ./your-path/min-working-example/

>>> pygna build-rwr-diffusion barabasi.interactome.tsv --output-file interactome_RWR.hdf5

>>> pygna test-topology-rwr barabasi.interactome.tsv disgenet_cancer_groups_subset.gmt interactome_RWR.hdf5 ../min_working_example/results/_topology_rwr.csv 
--number-of-permutations 1000 --cores 4

>>> pygna paint-datasets-stats interactome_table_RW.csv ./ example1

You can look at the plot of the results in the `example1_results.pdf` file, and the corresponding table in  `example1_table_RW.csv`.

Example 2: Running pygna GNA analysis
+++++++++++++++++++++++++++++++++++++

In the `min-working-example` we provide some basic files to run this minimal example. You can follow this instructions to run a network 
association analysis between two genesets.

>>> cd ./your-path/min-working-example/

skip this step if the matrix is already computed

>>> pygna build-rwr-diffusion barabasi.interactome.tsv --output-file interactome_RWR.hdf5

The association analysis is run N x M times (N number of genesets, M number of pathways), we use only 50 permutations in this example to avoid long computations; however, the recommended value is 1000.

>>> pygna test-association-rwr barabasi.interactome.tsv disgenet_cancer_groups_subset.gmt interactome_RWR.hdf5 table_association_rwr.csv -B disgenet_cancer_groups_subset.gmt --keep --number-of-permutations 100 --cores 4

If you don't include the --show-results flag at the comparison step, plot the matrix as follows

>>> pygna paint-comparison-matrix table_association_rwr.csv heatmap_association_rwr.png --rwr --annotate

( include the -rwr flag for the right color coding 
and --annotate for annotating the heatmap with the pvalue of each test ) 

The -k flag, keeps the -B geneset and permutes only on the set A.


WARNING: In this case, both A and B genesets are the same, usually you would use a different B geneset to check all the combinations.

If setname B is not passed, the analysis is run between each couple of setnames in the geneset.

>>> pygna test-association-rwr barabasi.interactome.tsv disgenet_cancer_groups_subset.gmt interactome_RWR.hdf5 table_within_comparison_rwr.csv --number-of-permutations 100 --cores 4 --show-results

>>> pygna paint-comparison-matrix table_within_comparison_rwr.csv heatmap_within_comparison_rwr.png --rwr --single-geneset



