import logging
import random
import networkx as nx
import numpy as np
import scipy
import pandas as pd
import sys
import matplotlib.pyplot as plt
import pygna.diagnostic as diag
import pygna.output as out
import pygna.parser as ps
import multiprocessing
import time
import seaborn as sns
from palettable.colorbrewer.diverging import *
from palettable.colorbrewer.sequential import *
import scipy.stats as stats


class Painter:
    def __init__(self, network, output_folder, name, diz={}):

        self.network = network
        self.output_folder = output_folder
        self.name = name
        self.diz = diz

        if (type(self.network) is nx.Graph) or (type(self.network) is nx.DiGraph):
            self.universe = set(self.network.nodes())
        elif type(self.network) is dict:
            self.universe = set(self.network.keys())
        else:
            logging.error("Unknown network type: %s" % type(self.network))
            sys.exit(-1)


class Painter_RW(Painter):
    def __init__(self, network, output_folder, name, RW_dict, geneset, diz={}):
        Painter.__init__(self, network, output_folder, name, diz={})
        self.RW_nodes = list(RW_dict["nodes"])
        self.geneset = list(geneset)
        self.geneset_index = [list(RW_dict["nodes"]).index(i) for i in geneset]
        self.RW_matrix = RW_dict["matrix"][self.geneset_index, :][:, self.geneset_index]

    def plot_matrix(self, show_labels=False):

        logging.info("Plotting figure as " + str(self.output_folder + self.name))

        fig, axes = plt.subplots(1, 1)
        fig.subplots_adjust(left=0.2, right=0.99, bottom=0.2, top=0.99)
        axes.set_title("Diffusion matrix of geneset's nodes")
        mask = np.eye(self.RW_matrix.shape[0])

        logging.info(
            str(np.array(self.geneset)[:, np.newaxis].shape) + str(self.RW_matrix.shape)
        )

        if show_labels:
            tab_conversion = pd.read_table(
                "/home/viola/Desktop/geneset-network-analysis/primary_data/entrez_name.tsv",
                sep="\t",
            )
            labels = []
            for i in self.geneset:
                name = tab_conversion[tab_conversion["EntrezID"] == int(i)][
                    "Symbol"
                ].values.tolist()
                if len(name) > 0:
                    labels.append(str(name[0]))
                else:
                    labels.append(i)

            g2 = sns.heatmap(
                self.RW_matrix - mask,
                cmap="OrRd",
                square=True,
                ax=axes,
                xticklabels=labels,
                yticklabels=labels,
                cbar=True,
                vmin=0,
                vmax=np.max(self.RW_matrix - mask),
            )
            g2.set_yticklabels(g2.get_yticklabels(), rotation=0, fontsize=6)
            g2.set_xticklabels(g2.get_xticklabels(), rotation=90, fontsize=6)
        else:
            g2 = sns.heatmap(self.RW_matrix, cmap="OrRd", mask=mask, ax=axes, cbar=True)

        fig.savefig(
            self.output_folder + "RW_matrix_" + self.name + ".pdf", format="pdf"
        )

    def restrict_network(self, geneset):

        self.__network = nx.DiGraph(
            self.__network.subgraph(max(nx.connected_components(self.__network), key=len))
        )

    def add_annotation_nodes_in_geneset(self, geneset, annotation="in_subset"):

        logging.info("adding annotation to nodes in geneset")
        logging.info("%d genes in geneset,%d mapped in network ",
            (len(geneset), len(self.__universe.intersection(set(geneset)))),
        )
        geneset = self.__universe.intersection(set(geneset))
        dict_nodes = {}
        for n in self.__network.nodes():
            if n in geneset[setname]:
                dict_nodes[n] = True
            else:
                dict_nodes[n] = False
        nx.set_node_attributes(self.__network, annotation, dict_nodes)

    def add_weights_edges(self, dict_edges):
        pass

    def draw_graphml(self):

        """
        Draws a .graphml file from the network
        """
        # mapping geneset
        logging.info("Drawing graph for " + str(name))
        nx.write_graphml(
            self.__network, self.__output_folder + self.__name + ".graphml"
        )


def paint_datasets_stats( table_filename: 'pygna results table',
                            output_file: 'figure file, use pdf or png extension',
                            alternative="greater"
    ):

    palette_binary = RdBu_4.mpl_colors[0::3]
    with open(table_filename, "r") as f:
        table = pd.read_csv(f, sep=",")

    stat_name = table["analysis"][0]
    n_permutations = table["number_of_permutations"][0]
    data = pd.DataFrame()
    if alternative == "greater":
        for i, row in table.sort_values(
            by=["empirical_pvalue", "observed"], ascending=[True, False]
        ).iterrows():
            data = data.append(
                {
                    "setname": row["setname"],
                    "pvalue": row["empirical_pvalue"],
                    "type": "observed",
                    "value": row["observed"],
                    "boundary_max": row["observed"] + 0,
                    "boundary_in": row["observed"] - 0,
                    "err": 0,
                },
                ignore_index=True,
            )
            data = data.append(
                {
                    "setname": row["setname"],
                    "type": "null",
                    "pvalue": row["empirical_pvalue"],
                    "value": row["mean(null)"],
                    "boundary_max": row["mean(null)"] + row["var(null)"],
                    "boundary_min": row["mean(null)"] - row["var(null)"],
                    "err": row["var(null)"],
                },
                ignore_index=True,
            )

        data = data.reset_index(drop=True)
        data["ind"] = data.index

        g = sns.catplot(
            x="value",
            y="setname",
            hue="type",
            data=data,
            height=6,
            kind="bar",
            palette=palette_binary,
            orient="h",
            aspect=1,
        )

        k = 0
        col = ["black", "#b64d50"]
        max_obs = np.max(data[data["type"] == "observed"]["value"])
        for i, r in data[
            data["type"] == "observed"
        ].iterrows():  # .sort_values(by=["pvalue"],ascending=True)
            if r["pvalue"] == 0:
                pval = "<%1.1E" % (1.0 / n_permutations)
            else:
                pval = "=%.4f" % r["pvalue"]
            g.facet_axis(0, 0).annotate(
                "pval:" + pval,
                xy=(r["value"], k),
                xytext=(r["value"], k),
                color=col[r["pvalue"] < (0.1 * 2 / len(data))],
                fontsize=6,
            )
            k += 1

        g.facet_axis(0, 0).xaxis.grid(color="gray", linestyle="dotted", linewidth=0.5)
        g.set_xlabels(stat_name)
        g.set_ylabels("")
        g.despine(bottom=True, left=True)

    else:

        for i, row in table.sort_values(
            by=["empirical_pvalue", "observed"], ascending=[True, True]
        ).iterrows():
            data = data.append(
                {
                    "setname": row["setname"],
                    "pvalue": row["empirical_pvalue"],
                    "type": "observed",
                    "value": row["observed"],
                    "boundary_max": row["observed"] + 0,
                    "boundary_in": row["observed"] - 0,
                    "err": 0,
                },
                ignore_index=True,
            )

            data = data.append(
                {
                    "setname": row["setname"],
                    "type": "null",
                    "pvalue": row["empirical_pvalue"],
                    "value": row["mean(null)"],
                    "boundary_max": row["mean(null)"] + row["var(null)"],
                    "boundary_min": row["mean(null)"] - row["var(null)"],
                    "err": row["var(null)"],
                },
                ignore_index=True,
            )

        # data=data.sort_values(by=["pvalue"],ascending=True)
        data = data.reset_index(drop=True)
        data["ind"] = data.index
        # data=data.iloc[0:20]

        # fig, ax = plt.subplots(1, 1, figsize=(10, 20), sharex=True,sharey=True)

        g = sns.catplot(
            x="value",
            y="setname",
            hue="type",
            data=data,
            height=6,
            kind="bar",
            palette=palette_binary,
            orient="h",
            aspect=1,
        )

        k = 0
        col = ["black", "#b64d50"]
        max_obs = np.max(data[data["type"] == "observed"]["value"])
        for i, r in data[
            data["type"] == "observed"
        ].iterrows():  # .sort_values(by=["pvalue"],ascending=True)
            if r["pvalue"] == 0:
                pval = "<%1.1E" % (1.0 / n_permutations)
            else:
                pval = "=%.4f" % r["pvalue"]
            g.facet_axis(0, 0).annotate(
                "pval:" + pval,
                xy=(r["value"], k),
                xytext=(r["value"] + 0.05, k - 0.10),
                color=col[r["pvalue"] < (0.1 * 2 / len(data))],
                fontsize=6,
            )
            k += 1

        g.facet_axis(0, 0).xaxis.grid(color="gray", linestyle="dotted", linewidth=0.5)
        g.set_xlabels(stat_name)
        g.set_ylabels("")
        g.despine(bottom=True, left=True)

    if output_file.endswith('.pdf'):
        plt.savefig(output_file, format="pdf")
    elif output_file.endswith('.png'):
        plt.savefig(output_file, format="png")
    else:
        logging.warning('The null distribution figure can only be saved in pdf or png, forced to png')
        plt.savefig(output_file+'.png', format="png")

def paint_comparison_matrix(table_filename: 'pygna comparison output', 
                            output_file: 'output figure file, specify png or pdf file',
                            rwr: 'use rwr is the table comes from a rwr analysis'= False,
                            single_geneset: 'use true if the comparison has been done for a single file'=False,
                            annotate: 'set true if uou want to print the pvalue inside the cell' = False):

    if rwr:
        palette = OrRd_9.mpl_colors
    else:
        palette = RdBu_11.mpl_colors
        
    with open(table_filename, "r") as f:
        table = pd.read_csv(f, sep=",")

    stat_name= table["analysis"][0]
    n_permutations = table["number_of_permutations"][0]

    if single_geneset:
        table=table.loc[:,['observed','setname_A','setname_B','empirical_pvalue']]
        for i in set(table['setname_A'].values.tolist()).union(set(table['setname_B'].values.tolist())):
            table=table.append({'setname_A':i, 'setname_B':i, 'observed':0, 'empirical_pvalue':1}, ignore_index=True)


    pivot_table = table.pivot(values="observed", index="setname_A", columns="setname_B")
    pivot_table = pivot_table.fillna(0)
    matrix = (pivot_table.T+pivot_table)-pivot_table.T*np.eye(len(pivot_table))
    if single_geneset:
        mask = np.triu(np.ones((len(pivot_table),len(pivot_table))))
    else:
        mask = np.ones((len(pivot_table),len(pivot_table)))-np.tril(np.ones((len(pivot_table),len(pivot_table))))

    if annotate:
        annot = table.pivot(
            values="empirical_pvalue", index="setname_A", columns="setname_B"
        )
        annot = annot.fillna(0)
        annot = (annot.T.values+annot.values)-annot.T.values*np.eye(len(annot))
    else:
        annot=False

    fig, axes = plt.subplots(1, 1, figsize=(10, 10))
    fig.subplots_adjust(left=0.3, right=0.99, bottom=0.3, top=0.99)
    
    if rwr:
        print('plotting rwr')
        g2 = sns.heatmap(
                matrix,
                cmap=palette,
                ax=axes,
                square=True,
                xticklabels=1,
                yticklabels=1,
                mask=mask,
                annot=annot,
                cbar=True,
                linewidths=0.1,
                linecolor="white",
            )  
    else:
        g2 = sns.heatmap(
                matrix,
                cmap=palette,
                ax=axes,
                square=True,
                xticklabels=1,
                yticklabels=1,
                mask=mask,
                annot=annot,
                center=0,
                cbar=True,
                linewidths=0.1,
                linecolor="white",
            )  
    g2.set_yticklabels(g2.get_yticklabels(), rotation=0, fontsize=6)
    g2.set_xticklabels(g2.get_xticklabels(), rotation=90, fontsize=6)
    axes.set_xlabel("")
    axes.set_ylabel("")

    if output_file.endswith('.pdf'):
        plt.savefig(output_file, format="pdf")
    elif output_file.endswith('.png'):
        plt.savefig(output_file, format="png")
    else:
        logging.warning('The null distribution figure can only be saved in pdf or png, forced to png')
        fig.savefig(output_file+'.png', format="png")


def plot_adjacency(
    network: "network_filename",
    output_file: 'use png or pdf for output figure', 
    clusters_file: "file of clusters to order the nodes" = None,
    size: "number of genes to plot, allows to cut the adjacency matrix" = None,
    ):

    """
    This function plots the adjacency matrix of a network. If a geneset file is passed,
    the matrix is organised by the different sets in the geneset.
    For the moment the genelist needs to be complete and non overlapping.
    """

    graph = ps.__load_network(network)
    if len(graph.nodes) > 1000:
        logging.warning("Graph is larger than 1k nodes, plotting might take too long")

    nodelist = None
    s = 0
    nodelabels = []
    if clusters_file:
        geneset = ps.__load_geneset(clusters_file)
        nodelist = [k for i, v in geneset.items() for k in v]
        for i, v in geneset.items():
            s += 1
            nodelabels.append([s] * len(v))
    nodelabels = [i for j in nodelabels for i in j]
    nodelabels = np.asarray(nodelabels)[np.newaxis, :]

    matrix = nx.adjacency_matrix(graph, nodelist=nodelist).toarray()

    if size:
        matrix = matrix[0 : int(size), 0 : int(size)]
        nodelabels = nodelabels[:, 0 : int(size)]

    if clusters_file:
        f, axes = plt.subplots(
            2,
            2,
            figsize=(10, 10),
            gridspec_kw={
                "width_ratios": [1, 100],
                "height_ratios": [1, 100],
                "wspace": 0.005,
                "hspace": 0.005,
            },
        )
        sns.heatmap(
            [[0, 0], [0, 0]],
            cmap="Greys",
            vmin=0,
            vmax=1,
            square=True,
            ax=axes[0, 0],
            cbar=False,
            xticklabels=False,
            yticklabels=False,
        )
        sns.heatmap(
            matrix,
            cmap="Greys",
            vmin=0,
            vmax=1,
            square=True,
            ax=axes[1, 1],
            cbar=False,
            xticklabels=False,
            yticklabels=False,
        )

        count = 0
        for i, v in geneset.items():
            if size and count < int(size):
                axes[0, 1].annotate(i, xy=(count, 0), xytext=(count, 0))
            elif size and count >= int(size):
                pass
            else:
                axes[0, 1].annotate(i, xy=(count, 0), xytext=(count, 0))
            count = count + len(v)

        g = sns.heatmap(
            nodelabels,
            xticklabels=False,
            yticklabels=False,
            vmin=1,
            vmax=len(geneset.keys()),
            square=False,
            cbar=False,
            linewidths=0,
            cmap="Set3",
            ax=axes[0, 1],
        )
        g = sns.heatmap(
            nodelabels.T,
            xticklabels=False,
            yticklabels=False,
            vmin=1,
            vmax=len(geneset.keys()),
            square=False,
            cbar=False,
            linewidths=0,
            cmap="Set3",
            ax=axes[1, 0],
        )
    else:
        f, axes = plt.subplots(1, 1, figsize=(10, 10))
        sns.heatmap(
            matrix,
            cmap="Greys",
            vmin=0,
            vmax=1,
            square=True,
            ax=axes[1, 1],
            cbar=False,
            xticklabels=False,
            yticklabels=False,
        )

    if output_file.endswith('.pdf'):
        plt.savefig(output_file, format="pdf")
    elif output_file.endswith('.png'):
        plt.savefig(output_file, format="png")
    else:
        logging.warning('The null distribution figure can only be saved in pdf or png, forced to png')
        f.savefig(output_file+'.png', format="png")