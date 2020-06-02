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

It is possible to implement custom functions that Pygna uses for the calculation of statistics.

.. autofunction:: statistical_test.StatisticalTest

Functions can be written anywhere but the return value of each function must be a `float`.

For example, it is possible to define a function, such as:

```def my_custom_function_test([...]) ->float:
[...]
return  value
```
In Pygna, it is possible to call the class constructor, passing the function as parameter:
`st_test = StatisticalTest(my_custom_function_test, **kwargs)`

Currently are implemented the following diffusion methods:

* geneset_localisation_statistic_median
* geneset_localisation_statistic
* geneset_module_statistic
* geneset_total_degree_statistic
* geneset_internal_degree_statistic
* geneset_RW_statistic

Stastistical Diffusion
++++++++++++++++++++++

As in the Statistical Test, also for the statistical diffusion class it is possible to define custom functions and use them in Pygna

.. autofunction:: statistical_diffusion.DiffusionTest

Functions can be written anywhere but the return value of each function must be a `float`.

For example, a function is defined as follows:

```
def my_custom_function_diffusion([...]) ->float:
[...]
return value
```

In Pygna, it is possible to call the class constructor, passing the function as parameter:
`t_test = DiffusionTest(my_custom_function_diffusion, **kwargs)`

Currently are implemented the following diffusion methods:

* weights_diffusion
* hotnet_diffusion


Statistical Comparison
++++++++++++++++++++++

The Statistical Comparison is another class where it is possible to use custom functions during the elaboration of Pygna.

.. autofunction:: statistical_comparison.StatisticalComparison

Functions can be written anywhere but the return value of each function must be a `float`.

For example, a function is defined as follows:

```
def my_custom_function_comparison([...]) ->float:
[...]
return value
```

In Pygna, it is possible to call the class constructor, passing the function as parameter:
`comp_test = StatisticalComparison(my_custom_function_comparison, **kwargs)`

Currently are implemented the following comparison methods:

* comparison_shortest_path
* comparison_random_walk

