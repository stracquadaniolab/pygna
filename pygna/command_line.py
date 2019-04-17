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
"""
autodoc
"""

logging.basicConfig(level=logging.INFO)

def main():
    """
    argh dispatch
    """
    argh.dispatch_commands([cmd.test_degree_distribution,
                            simulations.generate_simulated_network,
                            cmd.hierarchical_clustering_SP,
                            cmd.hierarchical_clustering_RW,
                            utils.new_network,
                            cmd.network_summary,
                            utils.convert_gmt,
                            paint.paint_datasets_stats,
                            paint.paint_comparison_stats,
                            paint.paint_comparison_RW,
                            cmd.analyse_internal_degree,
                            cmd.analyse_module,
                            cmd.analyse_location,
			                cmd.analyse_location_median,
                            cmd.analyse_RW,
                            cmd.comparison_shortest_path,
                            cmd.comparison_random_walk,
                            cmd.build_distance_matrix,
                            cmd.build_RWR_diffusion,
                            cmd.build_random_walk,
                            cmd.build_hitting_time_matrix])

if __name__ == "__main__":
    """
    MAIN
    """
    main()
