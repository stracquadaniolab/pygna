Statistical Diffusion
=====================

.. autoclass:: pygna.statistical_diffusion.DiffusionTest
    :members:

.. autofunction:: pygna.statistical_diffusion.map_table2matrix

.. autofunction:: pygna.statistical_diffusion.weights_diffusion_statistic

.. autofunction:: pygna.statistical_diffusion.hotnet_diffusion_statistic


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
Load the dictionary and the set the data you want to use
    >>> rw_dict = {"nodes": read_distance_matrix(rwr_matrix_filename, in_memory=in_memory)[0],
    ...            "matrix": read_distance_matrix(rwr_matrix_filename, in_memory=in_memory)[1]}
    >>> name_column = "gene_name"
    >>> weight_column = "stat"
    >>> table = pygna.reading_class.ReadCsv(geneset_file, column_to_fill=name_column).get_data()

Step 3
______
Define your custom function that implements your statistical test.
It must always return a `float` value.

For example here it is reported the calculation of hotnet 2 diffusion statistics.

.. code-block:: python

    def hotnet_diffusion_statistic(matrix, weights, geneset_index):
        weights = np.diagflat(weights.T)
        product = np.matmul(matrix, weights)
        prob = [product[i, j] for i in geneset_index for j in geneset_index if i != j]
        stat = np.sum(prob)
        return stat

Step 4
______
Pass your custom function to Pygna statistical class like:
    >>> st_test = pygna.statistical_diffusion.DiffusionTest(hotnet_diffusion_statistic, rw_dict["nodes"], rw_dict["matrix"], table, names_col=name_column, weights_col=weight_column)
