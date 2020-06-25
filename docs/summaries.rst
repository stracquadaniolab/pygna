Network Info
-----------------------------------

Before running any geneset-network analysis is a good practice to extract basic information on the network and the geneset or to visualise the network.
We provide a function to obtain the network summary ( or the summary of a single genesets) and another utility to write an annotated graphml file ( visualise it on Cytoscape ).

Network Properties
++++++++++++++++++++


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

.. autofunction:: pygna.command.network_graphml

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


    .. code-block:: text

    usage: pygna get-connected-components [-h] [-c] network-file geneset-file setname o graphml

    positional arguments:
      network-file          network file
      geneset-file          GMT geneset file
      setname               The setname to analyse
      o                     The output file name (should be gmt)
      graphml               The name of the graphml file

    optional arguments:
      -h, --help            show this help message and exit
      -c, --convert-entrez  Convert EntrezID->Symbol (default: True)

