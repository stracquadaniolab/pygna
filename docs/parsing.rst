File Parsing
============

The PyGNA package is able to read and write networks and datasetas in different formats.

Networks
--------

.. autofunction:: pygna.command.__load_network



Genesets
--------


.. autofunction:: pygna.command.__load_geneset


Large matrices
+++++++++++++

Along with the network and the geneset, some analyses require an additional large matrix to be passed as input.
In particular, the analysis of shortest path and diffusion are evaluating a matrix of shape NxN (N being the number of nodes),
since those are invariant to the geneset analysed, they must be evaluated and saved only once.

Shortest Paths matrix
---------------------

To evaluate the matrix of shortest paths we can use the following function:

.. autofunction:: pygna.command.build_distance_matrix

Diffusion matrix
---------------------

To evaluate the diffusion matrix there are different functions, depending on the type of
analysis we want to perform:

.. autofunction:: pygna.command.build_RWR_diffusion
.. autofunction:: pygna.command.build_random_walk
