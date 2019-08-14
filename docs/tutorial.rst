Tutorial
========


Complete Analysis of one geneset
--------------------------------

Shortest Paths matrix
+++++++++++++++++++++

To evaluate the matrix of shortest paths we can use the following function:


Diffusion matrix
++++++++++++++++

To evaluate the diffusion matrix we have the below function that
implements a Randowm Walk with Restart algorithm. The $beta$ parameter
is set to $0.80$ as default, but can be given by the user.



Pipelines
============

The package is integrable in Snakemake pipelines. We already provide
some sample analyses, but we encourage to explore all the functionalities
of the package.

Complete Analysis of one geneset
--------------------------------
a

Differential Expression Diffusion
---------------------------------
a

Pathway Analysis on Networks
----------------------------

a


Genesets creation and conversion
--------------------------------

The utils submodule provides a function to create the gmt files from .csv tables. In particular,
this function allows the user to filter the csv, for example, in case the table is the result of
a **differential expression** analysis, specifying the column and condition for the filter, the
user can obtain the gmt file of the differentially expressed terms.


Another utility allows to convert the gene names between gene-symbols and entrez codes. Usually, for consistency reasons,
genesets and networks are reported with entrez identifiers, while for interpretability scientists
prefer to read the names as gene symbols.

We provide a function to convert gmt names

whose core is the converter class

.. autoclass:: pygna.utils.Converter
