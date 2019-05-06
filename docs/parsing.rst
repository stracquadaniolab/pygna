Input files
============

The pyGNA package uses three main objects as input: genesets, gene tables and networks. 

Genesets
--------

Genesets are read in .gmt format through the following function

.. autofunction:: pygna.command.__load_geneset

The user can specify if all sets are read or can restrict the parser to return only the geneset with a specific setname.  

We provide functions for the creation of gmt files from tables and for the name conversion of gene_names.
Check the utilities section for further information.

Networks
--------

Networks are read in tsv format ( node_A \tab node_B ) through the function below

.. autofunction:: pygna.command.__load_network


Gene Tables
-----------

When performing the node weighting analysis, the user needs to specify not only the geneset
but the nodes for the entire network. In this case, the input can be a gene table, for 
example those returned by the deseq2 analysis. 

Large matrices
--------------

Along with the network and the geneset, some analyses require an additional large matrix to be passed as input.
In particular, the analysis of shortest path and diffusion are evaluating a matrix of shape NxN (N being the number of nodes),
since those are invariant to the geneset analysed, they must be evaluated and saved only once.

Shortest Paths matrix
+++++++++++++++++++++

To evaluate the matrix of shortest paths we can use the following function:

.. autofunction:: pygna.command.build_distance_matrix

Diffusion matrix
++++++++++++++++

To evaluate the diffusion matrix we have the below function that 
implements a Randowm Walk with Restart algorithm. The $beta$ parameter 
is set to $0.80$ as default, but can be given by the user.

.. autofunction:: pygna.command.build_RWR_diffusion
