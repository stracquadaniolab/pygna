API
===

**To be completed**

Input files
-----------

The pyGNA package uses three main objects as input: genesets, gene tables and networks.

Genesets
++++++++

Genesets are read in .gmt format through the following function

.. autofunction:: pygna.parser.__load_geneset

The user can specify if all sets are read or can restrict the parser to return only the geneset with a specific setname.

We provide functions for the creation of gmt files from tables and for the name conversion of gene_names.
Check the utilities section for further information.

Networks
+++++++++

Networks are read in tsv format ( node_A \tab node_B ) through the function below

.. autofunction:: pygna.parser.__load_network
