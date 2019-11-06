import matplotlib
matplotlib.use("TkAgg")

import logging
import argh
import sys
import numpy as np
import networkx as nx
import pygna.command as cmd
import pygna.painter as paint
import pygna.utils as utils
import pygna.KS_test as KS
import pygna.block_model as bm
import pygna.degree_model as dm
"""
autodoc
"""

logging.basicConfig(level=logging.INFO)

def main():
    """
    argh dispatch
    """
    argh.dispatch_commands([# network summary and graph file
                            cmd.network_summary,
                            cmd.network_graphml,
                            # geneset network topoly analyses
                            cmd.test_topology_total_degree,
                            cmd.test_topology_internal_degree,
                            cmd.test_topology_module,
                            cmd.test_topology_sp,
                            cmd.test_topology_rwr,
                            #cmd.test_degree_distribution,
                            cmd.test_diffusion_weights,
                            cmd.test_diffusion_hotnet,
                            # comparison analysis
                            cmd.test_association_sp,
                            cmd.test_association_rwr,
                            # building functions
                            cmd.build_distance_matrix,
                            cmd.build_rwr_diffusion,
                            # paint
                            paint.paint_datasets_stats,
                            paint.paint_comparison_matrix,
                            paint.plot_adjacency,
                            # utils
                            utils.convert_gmt,
                            utils.geneset_from_table,
                            utils.convert_csv_names,
                            utils.generate_group_gmt,
                             #simulations
                            dm.generate_hdn_network,
                            bm.generate_simulated_network,])

if __name__ == "__main__":
    """
    MAIN
    """
    main()
