Command Line Tool
======================

Network Summaries and Visualisation
-----------------------------------

Before running any geneset-network analysis is a good practice to extract basic information
on the network and the geneset or to visualise the network. We provide a function
to obtain the network summary ( or the summary of a single genesets) and another
utility to write an annotated graphml file ( visualise it on Cytoscape ).

.. autofunction:: pygna.command.network_summary

.. autofunction:: pygna.command.network_graphml


Geneset Network Topology Analysis
---------------------------------

Here are all the GNT analyses that pygna can perform.
For each of them a single geneset topology is tested with the specified test statistics.


GNT Module analysis
+++++++++++++++++++

.. autofunction:: pygna.command.test_topology_module

GNT Degree analysis
+++++++++++++++++++

.. autofunction:: pygna.command.test_topology_internal_degree
.. autofunction:: pygna.command.test_topology_total_degree

GNT Shortest Path Analysis
+++++++++++++++++++++++++++

.. autofunction:: pygna.command.test_topology_sp

GNT Random Walk Analysis
++++++++++++++++++++++++

.. autofunction:: pygna.command.test_topology_rwr

Geneset Network Association Analysis
------------------------------------

Here are all the GNA analyses that PyGNA can perform.
For each of them two sets are compared, one can pass:
- only one geneset: association between each set is computed
- two genesets: association between all the terms in the two genesets is computed

GNA Shortest Path
+++++++++++++++++++

.. autofunction:: pygna.command.test_association_sp

GNA Random Walk
+++++++++++++++++++

.. autofunction:: pygna.command.test_association_rwr


Weights Diffusion Analysis
--------------------------

**Please be aware that this function is still a beta version.**

.. autofunction:: pygna.command.test_diffusion_hotnet

Large Matrices creation
-----------------------

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

.. autofunction:: pygna.command.build_rwr_diffusion

Show Results
------------

.. autofunction:: pygna.painter.paint_datasets_stats
.. autofunction:: pygna.painter.paint_comparison_matrix

Utilities
---------

.. autofunction:: pygna.utils.convert_gmt
.. autofunction:: pygna.utils.geneset_from_table
.. autofunction:: pygna.utils.convert_csv
