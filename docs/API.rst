API
===

**To be completed**

Input files
-----------

The pyGNA package uses three main objects as input: genesets, gene tables and networks.

Genesets
++++++++

Genesets are read in .gmt format through the following function

.. autofunction:: rc.ReadGmt(geneset_file)

The user can specify if all sets are read or can restrict the parser to return only the geneset with a specific setname.

We provide functions for the creation of gmt files from tables and for the name conversion of gene_names.
Check the utilities section for further information.

Networks
+++++++++

Networks are read in tsv format ( node_A \tab node_B ) through the function below

.. autofunction:: pygna.parser.__load_network

Statistical Test
++++++++++++++++

 .. autoclass:: statistical_test.StatisticalTest class.

It is possible to use custom statistical function to be used into the StatisticalTest class.
Statistical functions can be written anywhere but the return value of each function must be a
.. code-block:: python
    float

For example, it is possible to define a function, such as:

.. code-block:: python
    def geneset_RW_statistic(network, geneset, diz={}, observed_flag=False):
        try:
            diz["matrix"]
        except KeyError:
            print("The dictionary doesnt have a matrix key")
            raise

        geneset_index = [diz["nodes"].index(i) for i in geneset]
        prob = [diz["matrix"][i, j] for i in geneset_index for j in geneset_index if i != j]
        prob = np.sum(prob)
        return prob

In Pygna, it is possible to call the class constructor, passing the function as parameter:

.. code-block:: python
    st_test = StatisticalTest(st.geneset_RW_statistic, network, rw_dict)

Currently are implemented the following diffusion methods:

* geneset_localisation_statistic_median
* geneset_localisation_statistic
* geneset_module_statistic
* geneset_total_degree_statistic
* geneset_internal_degree_statistic
* geneset_RW_statistic

Statistical Diffusion
++++++++++++++++++++++

.. autoclass:: statistical_diffusion.DiffusionTest

It is possible to use custom statistical diffusion functions to be used into the StatisticalDiffusion class.
Statistical functions can be written anywhere but the return value of each function must be a
.. code-block:: python
    float

For example, a function is defined as follows:

.. code-block:: python
    def hotnet_diffusion_statistic(matrix, weights, geneset_index, diz={}, observed_flag=False):
        weights = np.diagflat(weights.T)
        if matrix.shape[1] != weights.shape[0]:
            logging.warning("pass the right shape for the weights")
        try:
            product = np.matmul(matrix, weights)
        except:
            logging.warning("error in matrix multiplication")

        prob = [product[i, j] for i in geneset_index for j in geneset_index if i != j]
        stat = np.sum(prob)
        return stat

It is possible to call the class constructor as:
.. code-block:: python
    t_test = DiffusionTest(sd.hotnet_diffusion_statistic, rw_dict["nodes"], rw_dict["matrix"], table, names_col=name_column, weights_col=weight_column)

Currently are implemented the following diffusion methods:

* weights_diffusion
* hotnet_diffusion


Statistical Comparison
++++++++++++++++++++++

.. autoclass:: statistical_comparison.StatisticalComparison

It is possible to use custom statistical comparison functions to be used into the StatisticalComparison class.
Functions can be written anywhere but the return value of each function must be a
.. code-block:: python
    float

For example, a function is defined as follows:

.. code-block:: python
    def comparison_shortest_path(network, genesetA, genesetB, diz):
        n = np.array([diz["nodes"].index(i) for i in genesetA])
        m = np.array([diz["nodes"].index(i) for i in genesetB])
        if len(n) == 0 or len(m) == 0:
            logging.info("Geneset length is equal to 0")
            sys.exit()
        cum_sum = calculate_sum(n, m, diz)
        cum_sum += calculate_sum(m, n, diz)

        d_AB = cum_sum / float(len(genesetA) + len(genesetB))

        d_A = st.geneset_localisation_statistic(network, genesetA, diz)
        d_B = st.geneset_localisation_statistic(network, genesetB, diz)
        return d_AB - (d_A + d_B) / 2


It is possible to call the class constructor as:

.. code-block:: python
    comp_test = StatisticalComparison(sc.comparison_shortest_path, network, diz=sp_diz, n_proc=cores)

Currently are implemented the following comparison methods:

* comparison_shortest_path
* comparison_random_walk

