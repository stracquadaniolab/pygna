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


**#TODO: add a usage example test_centrality.py**
For the new centrality test create a single py file with a


.. toctree::
    :maxdepth: 1

    vignettes



Diagnostic
+++++++++++++++++++

For the GNT tests we also provide some diagnostic tools.

**#TODO: Add example of distribution plot**
