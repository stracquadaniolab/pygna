Simulated Networks
==================

We have created a framework for testing the performance of Geneset Network analysis: it allows
the experimenter to create simulated networks specifying many parameters.
The model we use is the Stochastic Block Model and Degree Model that provides a way to create the network
and identify clusters into it.

The function to generate a simulated dataset is:

.. autofunction:: pygna.degree_model.generate_hdn_network

.. autofunction:: pygna.block_model.generate_sbm_network

.. autofunction:: pygna.painter.plot_adjacency

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
-------------------------------------------------

**Generation of the network and genesets**

>>> pygna generate-simulated-network ../simulationBM.yaml

