Benchmarking
==================

We have created a framework for testing the performance of Geneset Network
analysis: it allows the experimenter to create simulated networks specifying
many parameters. The model we use is the Stochastic Block Model and Degree Model
that provides a way to create the network and identify clusters into it.

Stochastic Block Model
--------------------------

We provide both functions to generate stochastic block model simualated data like below.
On the left is the model used for GNT benchmarking and on the right the one for GNA benchmarking.

For a quick primer on SBM, you can check the page below:

.. toctree::
    :maxdepth: 1

    show_sbm

For each simulation we return both the whole network and a gmt file with the blocks.

.. image:: _static/figure_sim_examples.png
   :width: 600

GNT Benchmarking
+++++++++++++++++++

.. autofunction:: pygna.block_model.generate_gnt_sbm


GNA Benchmarking
+++++++++++++++++++

.. autofunction:: pygna.block_model.generate_gna_sbm




High Degree Nodes simulations
-------------------------------


Generate the network and geneset file
++++++++++++++++++++++++++++++++++++++++


The high degree nodes (HDN) model generates networks with a controllable number
of hubs, HDNs, whose probability of connection with another node
is higher than the baseline probability assigned to any other node in the
network.

.. autofunction:: pygna.degree_model.generate_hdn_network


Add genesets
++++++++++++++++++++++++++++++++++++++++

Given the generated network and node list of HDNs, we can then generate
novel genesets made of mixtures of the two.

We show the idea of how they are generated below. First is the
original network with a number of HDNs, then the partial, extended, and branching genesets.

.. image:: _static/Simulations_HDN_concept.png
   :width: 600

Add extended genesets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: pygna.degree_model.hdn_add_extended


Add partial genesets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Partial genesets are

.. autofunction:: pygna.degree_model.hdn_add_partial


Add branching genesets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: pygna.degree_model.hdn_add_branching





General model (old)
++++++++++++++++++++

This function is still working, however now networkx provides a function
to generate SBM networks. We will then change this function to use
directly networkx for the generation.

All the parameters can be specified in a yaml file that is passed as input. As output,
we obtain a file with the network and a .gmt file where the nodes have been grouped
by the respective cluster they are in.


Example yaml
------------

.. code-block:: yaml

    BlockModel :
      n_nodes: 100
      matrix: [[0.8,0.1,0.1],[0.1,0.8,0.1],[0.1,0.1,0.8]]
      nodes: ""
      nodes_in_cluster: None
    Simulations:
      n_simulated: 1
      output_folder: "/home/viola/Desktop/geneset-network-analysis/processed_data/simulated_data/"
      suffix: attempt1


The **Simulations** parameters are used to specify the number and name of the output:

+ ``n_simulated`` we can specify the number of networks we want to generate with the same settings
+ ``output_folder`` specifies the folder where the output files are going to be saved
+ ``suffix`` is the suffix used for the output.

For each simulated datasets there are two output files:
- suffix_network_$simulation.tsv : the network file
- suffix_genes_$simulation.gmt : the gene list grouped by cluster


The **BlockModel** parameters are those used to generate the SBM:

+  ``n_nodes`` number of nodes of the network
+  ``matrix`` SBM matrix
+  ``nodes: ""`` names of the nodes, if not specified N$number
+  ``nodes_in_cluster: None``: Number of nodes thatwe want to assign to each cluster


Example simulated dataset generation and Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Generation of the network and genesets**

.. code-block:: bash

    $ pygna generate-simulated-network ../simulationBM.yaml
