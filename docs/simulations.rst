Simulated Networks
==================

We have created a framework for testing the performance of Geneset Network analysis: it allows
the experimenter to create simulated networks specifying many parameters.
The model we use is the Stochastic Block Model that provides a way to create the network
and identify clusters into it.

The function to generate a simulated dataset is:

.. autofunction:: pygna.simulations.generate_simulated_network


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

**Generation of the network and SP/RWR matrices**

.. code-block:: bash

  pygna generate-simulated-network "/home/viola/Desktop/geneset-network-analysis/analyses/simulationBM.yaml"

  pygna build-distance-matrix '/home/viola/Desktop/geneset-network-analysis/processed_data/simulated_data/attempt1_network_0.tsv' '/home/viola/Desktop/geneset-network-analysis/processed_data/simulated_data/attempt1_network_0_SP.hdf5'

  pygna build-RWR-diffusion '/home/viola/Desktop/geneset-network-analysis/processed_data/simulated_data/attempt1_network_0.tsv' --output-file '/home/viola/Desktop/geneset-network-analysis/processed_data/simulated_data/attempt1_network_0_RWR'

**Analysis**

.. code-block:: bash

  #!/bin/bash

  NETWORK='/home/viola/Desktop/geneset-network-analysis/processed_data/simulated_data/attempt1_network_0.tsv'
  GMTFILE='/home/viola/Desktop/geneset-network-analysis/processed_data/simulated_data/attempt1_genes_0.gmt'
  OUTPATH='/home/viola/Desktop/geneset-network-analysis/processed_data/simulated_data/'
  DISTANCEM='/home/viola/Desktop/geneset-network-analysis/processed_data/simulated_data/attempt1_network_0_SP.hdf5'
  DIFFUSIONM='/home/viola/Desktop/geneset-network-analysis/processed_data/simulated_data/attempt1_network_0_RWR.pickle'
  NP=1000
  CORES=1

  pygna analyse-module $NETWORK  $GMTFILE $OUTPATH --number-of-permutations $NP --cores $CORES --create-output-LCC --show-results
  pygna analyse-location $NETWORK $DISTANCEM $GMTFILE $OUTPATH --number-of-permutations $NP --cores $CORES --show-results
  pygna analyse-RW $NETWORK $GMTFILE $DIFFUSIONM $OUTPATH --number-of-permutations $NP --cores $CORES --show-results
