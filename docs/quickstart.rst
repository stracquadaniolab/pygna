Quickstart
============

Installation
------------

The easiest and fastest way to install `pygna` using `conda`:

.. code-block:: bash

    $ conda install -c stracquadaniolab -c bioconda -c conda-forge pygna

Alternatively you can install it through `pip`:

.. code-block:: bash

    $ pip install pygna

We have prepared also a docker image that can be found at the following link:

https://github.com/stracquadaniolab/pygna/packages

By downloading the docker image named "Pygna", you will download a virtual machine with the latest version of pygna and its requirements.
First download and install the docker environment. Then run it either by command line or by opening the application.
We have setup the docker image in order to let you use it "out of the box" as stand-alone app, so just by running the following line:

.. code-block:: bash

    $ docker run docker.pkg.github.com/stracquadaniolab/pygna/pygna:latest

you will be prompt to the Pygna help section.

You can use pygna commands I/O functions using the classic docker functionality.
Here below is an example that paint the dataset statistics:

.. code-block:: bash

    $ docker run --rm -v "$PWD:$PWD" -w "$PWD" docker.pkg.github.com/stracquadaniolab/pygna/pygna:latest paint-datasets-stats pygna_gnt_results.csv pygna_gnt.png


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

**WARNING**: Pay attention to the fact we set the number of permutations to 1000, if you are just willing to test the behaviour, use 50 or 100
to speed up the process*

.. code-block:: bash

    $ cd ./your-path/min-working-example/

    $ pygna build-rwr-diffusion barabasi.interactome.tsv --output-file interactome_RWR.hdf5

    $ pygna test-topology-rwr barabasi.interactome.tsv disgenet_cancer_groups_subset.gmt interactome_RWR.hdf5 table_topology_rwr.csv --number-of-permutations 1000 --cores 4

    $ pygna paint-datasets-stats table_topology_rwr.csv barplot_rwr.png

You can look at the plot of the results in the `barplot_rwr.png` file, and the corresponding table in  `table_topology_rwr.csv`.

Example 2: Running pygna GNA analysis
+++++++++++++++++++++++++++++++++++++

In the `min-working-example` we provide some basic files to run this minimal example. You can follow this instructions to run a network
association analysis between two genesets.

.. code-block:: bash

    $ cd ./your-path/min-working-example/

Skip this step if the matrix is already computed

.. code-block:: bash

    $ pygna build-rwr-diffusion barabasi.interactome.tsv --output-file interactome_RWR.hdf5

The association analysis is run N x M times (N number of genesets, M number of pathways), we use only 50 permutations in this example to avoid long computations; however, the recommended value is 1000.

.. code-block:: bash

    $ pygna test-association-rwr barabasi.interactome.tsv disgenet_cancer_groups_subset.gmt interactome_RWR.hdf5 table_association_rwr.csv -B disgenet_cancer_groups_subset.gmt --keep --number-of-permutations 100 --cores 4

If you don't include the --results-figure flag at the comparison step, plot the matrix as follows

.. code-block:: bash

    $ pygna paint-comparison-matrix table_association_rwr.csv heatmap_association_rwr.png --rwr --annotate

( include the -rwr flag for the right color coding
and --annotate for annotating the heatmap with the pvalue of each test )

The -k flag, keeps the -B geneset and permutes only on the set A.


**WARNING**: In this case, both A and B genesets are the same, usually you would use a different B geneset to check all the combinations.

If setname B is not passed, the analysis is run between each couple of setnames in the geneset.

.. code-block:: bash

    $ pygna test-association-rwr barabasi.interactome.tsv disgenet_cancer_groups_subset.gmt interactome_RWR.hdf5 table_within_comparison_rwr.csv --number-of-permutations 100 --cores 4
    $ pygna paint-comparison-matrix table_within_comparison_rwr.csv heatmap_within_comparison_rwr.png --rwr --single-geneset



