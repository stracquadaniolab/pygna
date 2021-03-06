
Large Matrices creation
-----------------------

Along with the network and the geneset, some analyses require an additional large matrix to be passed as input.
In particular, the analysis of shortest path and diffusion are evaluating a matrix of shape :math:`NxN` (N being the number of nodes),
since those are invariant to the geneset analysed, they must be evaluated and saved only once.


Shortest Paths matrix
+++++++++++++++++++++


Build a shortest path distance matrix for a given network. Matrix can be saved as a .txt file or a .hdf5 one.

.. code-block:: text

    usage: pygna build-distance-matrix [-h] [-g] network-file output-file

    positional arguments:
      network-file          network file
      output-file           distance matrix output file, use .hdf5

    optional arguments:
      -h, --help            show this help message and exit
      -g, --giant-component-only
                            compute the shortest paths only for nodes in the giant component (default: True)


Diffusion matrix
++++++++++++++++

To evaluate the diffusion matrix we have the below function that implements a Randowm Walk with Restart algorithm.
The :math:`beta` parameter is set to 0.80 as default, but can be given by the user.


.. code-block:: text

    usage: pygna build-rwr-diffusion [-h] [-b BETA] [-o OUTPUT_FILE] network-file

    positional arguments:
      network-file          network file

    optional arguments:
      -h, --help            show this help message and exit
      -b BETA, --beta BETA  0.85
      -o OUTPUT_FILE, --output-file OUTPUT_FILE
                            distance matrix output file (use .hdf5) (default: -)



Converting tables and names
----------------------------------

Dataset from table
+++++++++++++++++++++

Converts a csv file to a GMT allowing to filter the elements using the values of one of the columns.
The user can specify the column used to retrieve the name of the objects and the filter condition.
The output can be either a GMT with the names of the genes that pass the filter or a csv with the whole filtered table, otherwise both can be created.

.. code-block:: text

    usage: pygna geneset-from-table [-h] [--output-gmt OUTPUT_GMT] [--output-csv OUTPUT_CSV] [-n NAME_COLUMN] [-f FILTER_COLUMN] [-a ALTERNATIVE] [-t THRESHOLD]
                                [-d DESCRIPTOR]
                                input-file setname

    positional arguments:
      input-file            input csv file
      setname               name of the set

    optional arguments:
      -h, --help            show this help message and exit
      --output-gmt OUTPUT_GMT
                            output gmt name (default: -)
      --output-csv OUTPUT_CSV
                            output csv name (default: -)
      -n NAME_COLUMN, --name-column NAME_COLUMN
                            column with the names (default: 'Unnamed: 0')
      -f FILTER_COLUMN, --filter-column FILTER_COLUMN
                            column with the values to be filtered (default: 'padj')
      -a ALTERNATIVE, --alternative ALTERNATIVE
                            alternative to use for the filter, with less the filter is applied <threshold, otherwise >= threshold (default: 'less')
      -t THRESHOLD, --threshold THRESHOLD
                            threshold for the filter (default: 0.01)
      -d DESCRIPTOR, --descriptor DESCRIPTOR
                            descriptor for the gmt file (default: -)


Convert gene names
+++++++++++++++++++++

`convert-gmt` is used to convert a GMT file, adding information about the Entrez ID or the symbol

.. code-block:: text

    usage: pygna convert-gmt [-h] [-e ENTREZ_COL] [-s SYMBOL_COL] gmt-file output-gmt-file conversion converter-map-filename

    positional arguments:
      gmt-file              gmt file to be converted
      output-gmt-file       output file
      conversion            e2s or s2e
      converter-map-filename
                            tsv table used to convert gene names

    optional arguments:
      -h, --help            show this help message and exit
      -e ENTREZ_COL, --entrez-col ENTREZ_COL
                            name of the entrez column (default: 'NCBI Gene ID')
      -s SYMBOL_COL, --symbol-col SYMBOL_COL
                            name of the symbol column (default: 'Approved symbol')


`generate-group-gmt` generates a GMT file of multiple setnames.
From the table file, it groups the names in the group_col (the column you want to use to group them) and prints the genes in the name_col.
Set the descriptor according to your needs

.. code-block:: text

    usage: pygna generate-group-gmt [-h] [-n NAME_COL] [-g GROUP_COL] [-d DESCRIPTOR] input-table output-gmt

    positional arguments:
      input-table           table to get the geneset from
      output-gmt            output GMT file

    optional arguments:
      -h, --help            show this help message and exit
      -n NAME_COL, --name-col NAME_COL
                            'Gene'
      -g GROUP_COL, --group-col GROUP_COL
                            'Cancer'
      -d DESCRIPTOR, --descriptor DESCRIPTOR
                            'cancer_genes'

`convert-csv`  is used to add a column with the entrezID or Symbols to a CSV file.

.. code-block:: text

    usage: pygna convert-csv [-h] [--converter-map-filename CONVERTER_MAP_FILENAME] [--output-file OUTPUT_FILE] [-e ENTREZ_COL] [-s SYMBOL_COL]
                         csv-file conversion original-name-col new-name-col geneset

    positional arguments:
      csv-file              csv file where to add a name column
      conversion            e2s or s2e
      original-name-col     column name to be converted
      new-name-col          name of the new column with the converted names
      geneset               the geneset to convert

    optional arguments:
      -h, --help            show this help message and exit
      --converter-map-filename CONVERTER_MAP_FILENAME
                            tsv table used to convert gene names (default: 'entrez_name.tsv')
      --output-file OUTPUT_FILE
                            if none, table is saved in the same input file (default: -)
      -e ENTREZ_COL, --entrez-col ENTREZ_COL
                            name of the entrez column (default: 'NCBI Gene ID')
      -s SYMBOL_COL, --symbol-col SYMBOL_COL
                            name of the symbol column (default: 'Approved symbol')


