.. _visualisation:

Visualisation
----------------

Pygna offers a range of plots to visualise the results of the analysis.

GNT plots
++++++++++++

GNT barplot
^^^^^^^^^^^^^^


**#TODO: add one plot and brief explaination**


.. autofunction:: pygna.painter.paint_datasets_stats

.. code-block:: text

    usage: pygna paint-datasets-stats [-h] [-a ALTERNATIVE] table-filename output-file

    positional arguments:
      table-filename        pygna results table
      output-file           figure file, use pdf or png extension

    optional arguments:
      -h, --help            show this help message and exit
      -a ALTERNATIVE, --alternative ALTERNATIVE
                            'greater'

GNT summary
^^^^^^^^^^^^^^


**#TODO: add one summary and brief explaination**

.. autofunction:: pygna.painter.paint_summary_gnt

.. code-block:: text

    usage: pygna paint-summary-gnt [-h] [-s SETNAME] [-t THRESHOLD] [-c COLUMN_FILTER] [--larger] [--less-tests LESS_TESTS] output-figure [input_tables [input_tables ...]]

    positional arguments:
      output-figure         output figure filename
      input_tables          -

    optional arguments:
      -h, --help            show this help message and exit
      -s SETNAME, --setname SETNAME
                            name of the dataset (default: -)
      -t THRESHOLD, --threshold THRESHOLD
                            Value to threshold the colors (default: 0.05)
      -c COLUMN_FILTER, --column-filter COLUMN_FILTER
                            column where the threshold is applied (default: 'empirical_pvalue')
      --larger              if True the threshold is set as lower limit (default: False)
      --less-tests LESS_TESTS
                            comma separated string of the tests that are significant if lower than expected, otherwise pass empty string (default: 'topology_sp')



GNA plots
++++++++++++

GNA heatmap
^^^^^^^^^^^^^^^^^^^^^^^

**#TODO: put two examples one for triangular heatmap comparison, one for whole heatmap association**

.. autofunction:: pygna.painter.paint_comparison_matrix

.. code-block:: text

    usage: pygna paint-comparison-matrix [-h] [-r] [-s] [-a] table-filename output-file

    positional arguments:
      table-filename        pygna comparison output
      output-file           output figure file, specify png or pdf file

    optional arguments:
      -h, --help            show this help message and exit
      -r, --rwr             use rwr is the table comes from a rwr analysis (default: False)
      -s, --single-geneset  use true if the comparison has been done for a single file (default: False)
      -a, --annotate        set true if uou want to print the pvalue inside the cell (default: False)






GNA association volcano
^^^^^^^^^^^^^^^^^^^^^^^

**#TODO: add volcano rwr and sp**

.. autofunction:: pygna.painter.paint_volcano_plot

.. code-block:: text

    usage: pygna paint-volcano-plot [-h] [-r] [-i ID_COL] [--threshold-x THRESHOLD_X] [--threshold-y THRESHOLD_Y] [-a] table-filename output-file

    positional arguments:
      table-filename        pygna comparison output
      output-file           output figure file, specify png or pdf file

    optional arguments:
      -h, --help            show this help message and exit
      -r, --rwr             use rwr is the table comes from a rwr analysis (default: False)
      -i ID_COL, --id-col ID_COL
                            'setname_B'
      --threshold-x THRESHOLD_X
                            0
      --threshold-y THRESHOLD_Y
                            2
      -a, --annotate        False

