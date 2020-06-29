Network Info
-----------------------------------

Before running any geneset-network analysis is a good practice to extract basic information on the network and the geneset or to visualise the network.
We provide a function to obtain the network summary ( or the summary of a single genesets) and another utility to write an annotated graphml file ( visualise it on Cytoscape ).

Network Properties
++++++++++++++++++++

`network-summary` saves the principal info of a graph:
- network properties
- degree distribution
- connected components diagnostic

If a geneset/setname is passed to the function, the properties of the subgraph are evaluated

.. code-block:: text

    pygna network-summary [-h] [-g GENESET_INPUT_FILE] [-s SETNAME] network-file text-output degree-figure-file c-components-figure-file

    positional arguments:
      network-file          network file
      text-output           output text file for the summary
      degree-figure-file    pdf or png file for the degree distribution
      c-components-figure-file
                            pdf or png file for the connected components distribution

    optional arguments:
      -h, --help            show this help message and exit
      -g GENESET_INPUT_FILE, --geneset-input-file GENESET_INPUT_FILE
                            geneset file (default: -)
      -s SETNAME, --setname SETNAME
                            specify a single geneset (default: -)


Cytoscape visualisation
++++++++++++++++++++++++

`network-graphml` generates a graphml file with nodes annotation.
Given a geneset, with k setnames, each node has k False/True annotations for each set.

Warning: without minimal, this function saves the full network.
The minimal graph saves only the nodes in the geneset and those that connect them with a shortest path.

.. code-block:: text

    usage: pygna network-graphml [-h] [-s SETNAME] [--giant-component-only] [-m] network-file geneset-file output-file

    positional arguments:
      network-file          network file
      geneset-file          geneset file
      output-file           graphml file for network for visualisation

    optional arguments:
      -h, --help            show this help message and exit
      -s SETNAME, --setname SETNAME
                            set name (default: -)
      --giant-component-only
                            saves only the giant component of the network (default: True)
      -m, --minimal         saves only the minimal graph (default: False)

Connected components
+++++++++++++++++++++

`get-connected-components` evaluate all the connected components in the subgraph pf the network with a given setname.
Multiple setnames can be passed to this function to analyze all of them in a run.
The file produces a GMT output and optionally a plot of the subnetwork with the connected components analysed.

Please notice that to convert the entrezID into Symbols, a stable internet connection is required

.. code-block:: text

        usage: pygna get-connected-components [-h] [--geneset-file GENESET_FILE]
                                            [-s SETNAME] [--graphml GRAPHML]
                                            [-t THRESHOLD] [-c]
                                            network-file output-gmt name

            This function evaluate all the connected components in the subgraph pf the network with a given setname.
            Multiple setnames can be passed to this function to analyze all of them in a run.
            The file produces a GMT output and optionally a plot of the subnetwork with the connected components analysed.
            Please notice that to convert the entrezID into Symbols, a stable internet connection is required


        positional arguments:
        network-file          network tsv file
        output-gmt            The output file name (should be gmt)
        name                  pass a name for the putput gmt terms

        optional arguments:
        -h, --help            show this help message and exit
        --geneset-file GENESET_FILE
                                GMT of the geneset file, is a file is passed please
                                add the setname (default: -)
        -s SETNAME, --setname SETNAME
                                The setname to analyse (default: -)
        --graphml GRAPHML     Pass a graphml filename to show the results on
                                Cytoscape (default: -)
        -t THRESHOLD, --threshold THRESHOLD
                                ignores all CC smaller than this value (default: 1)
        -c, --convert-entrez  pass flag to convert EntrezID->Symbol (default: False)

