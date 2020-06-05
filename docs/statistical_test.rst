Statistical Test
================

.. autoclass:: pygna.statistical_test.StatisticalTest
    :members:

.. autofunction:: pygna.statistical_test.geneset_localisation_statistic_median

.. autofunction:: pygna.statistical_test.geneset_localisation_statistic

.. autofunction:: pygna.statistical_test.geneset_module_statistic

.. autofunction:: pygna.statistical_test.geneset_total_degree_statistic

.. autofunction:: pygna.statistical_test.geneset_internal_degree_statistic

.. autofunction:: pygna.statistical_test.geneset_RW_statistic


Customizing the statistics
++++++++++++++++++++++++++

In addition to the provided functions, the user can define and use any function.
Follow the example below.

Step 1
______
Import Pygna and its libraries
    >>> from pygna import *

Step 2
______
Load the network you want to use
    >>> network = pygna.reading_class.ReadTsv("mynetwork.tsv").get_network()

Step 3
______
Define your custom function that implements your statistical test.
It must always return a `float` value.

For example here it is reported the calculation of the geneset total degree

.. code-block:: python

    def geneset_total_degree_statistic(network, geneset):
        degree = nx.degree(network)
        geneset = list(geneset)
        total = np.array([degree[g] for g in geneset])
        return np.average(total)

Step 4
______
Pass your custom function to Pygna statistical class like:
    >>> st_test = pygna.statistical_test.StatisticalTest(geneset_total_degree_statistic, network)

