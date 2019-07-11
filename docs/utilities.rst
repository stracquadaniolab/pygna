Utilities
============

Along with the main functions, this package provides some utilities to help with the practical usage of
pyGNA. 

Network Summary
---------------

Function to create a summary of the network properties. 

.. autofunction:: pygna.command.network_summary


Genesets creation and conversion
--------------------------------

The utils submodule provides a function to create the gmt files from .csv tables. In particular, 
this function allows the user to filter the csv, for example, in case the table is the result of 
a **differential expression** analysis, specifying the column and condition for the filter, the
user can obtain the gmt file of the differentially expressed terms.

.. autofunction:: pygna.utils.geneset_from_table

Another utility allows to convert the gene names between gene-symbols and entrez codes. Usually, for consistency reasons, 
genesets and networks are reported with entrez identifiers, while for interpretability scientists
prefer to read the names as gene symbols.

We provide a function to convert gmt names

.. autofunction:: pygna.utils.convert_gmt 

whose core is the converter class 

.. autoclass:: pygna.utils.Converter

Drawing Graphml files for Cytoscape
-----------------------------------


