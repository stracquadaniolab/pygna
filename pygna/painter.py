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
        print(np.min(self.RW_matrix))

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
            self.__network.subgraph(max(nx.connected_components(network), key=len))
        )

    def add_annotation_nodes_in_geneset(self, geneset, annotation="in_subset"):

        logging.info("adding annotation to nodes in geneset")
        logging.info(
            "%d genes in geneset,%d mapped in network ",
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


def paint_datasets_stats(
    table_filename, output_folder, stat_name, alternative="greater"
):

    palette_binary = RdBu_4.mpl_colors[0::3]
    with open(table_filename, "r") as f:
        table = pd.read_table(f, sep=",")
    print(table.head())
    # table=table.sort_values(by=["empirical_pvalue","observed"],ascending=[True, False])

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

        # Example data
        print(data["err"].values[:, np.newaxis].shape)

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
        # Example data
        print(data["err"].values[:, np.newaxis].shape)

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

    plt.savefig(output_folder + stat_name + "_results.pdf", format="pdf")
    # plt.savefig(output_folder+"PA_fisher_"+dataset_name+".jpg", format='jpg')


def paint_comparison_stats(table_filename, output_folder, stat_name):

    palette = RdBu_11.mpl_colors  # [0::3]
    with open(table_filename, "r") as f:
        table = pd.read_table(f, sep=",")
    print(table.head())

    n_permutations = table["number_of_permutations"][0]

    pivot_table = table.pivot(values="observed", index="setname_A", columns="setname_B")

    annot = table.pivot(
        values="empirical_pvalue", index="setname_A", columns="setname_B"
    ).values

    fig, axes = plt.subplots(1, 1, figsize=(10, 10))
    fig.subplots_adjust(left=0.3, right=0.99, bottom=0.3, top=0.99)

    pivot_table = pivot_table.fillna(0)

    g2 = sns.heatmap(
        pivot_table.T,
        cmap=palette,
        ax=axes,
        square=True,
        xticklabels=1,
        yticklabels=1,
        cbar=True,
        center=0,
        linewidths=0.1,
        linecolor="white",
    )  # ,annot=annot,fmt=""

    g2.set_yticklabels(g2.get_yticklabels(), rotation=0, fontsize=6)
    g2.set_xticklabels(g2.get_xticklabels(), rotation=90, fontsize=6)
    axes.set_xlabel("")
    axes.set_ylabel("")

    fig.savefig(output_folder + stat_name + "_comparison_heatmap.pdf", format="pdf")
    # fig.savefig(output_folder+"fig_comparison_heatmap.jpeg", format='jpeg')


def paint_comparison_RW(table_filename, output_folder, stat_name, single_geneset=False):

    palette = OrRd_9.mpl_colors  # [0::3]
    with open(table_filename, "r") as f:
        table = pd.read_table(f, sep=",")
    print(table.head())

    n_permutations = table["number_of_permutations"][0]

    if single_geneset:
        table=table.loc[:,['observed','setname_A','setname_B','empirical_pvalue']]
        for i in set(table['setname_A'].values.tolist()).union(set(table['setname_B'].values.tolist())):
            table=table.append({'setname_A':i, 'setname_B':i, 'observed':0, 'empirical_pvalue':1}, ignore_index=True)
        print(table)

    pivot_table = table.pivot(values="observed", index="setname_A", columns="setname_B")

    annot = table.pivot(
        values="empirical_pvalue", index="setname_A", columns="setname_B"
    ).values
    if single_geneset:

        fig, axes = plt.subplots(1, 1, figsize=(10, 10))
        fig.subplots_adjust(left=0.3, right=0.99, bottom=0.3, top=0.99)

        pivot_table = pivot_table.fillna(0)
        mask = np.triu(np.ones((len(pivot_table),len(pivot_table))))

        g2 = sns.heatmap(
            (pivot_table.T+pivot_table)/2,
            cmap=palette,
            ax=axes,
            square=True,
            xticklabels=1,
            yticklabels=1,
            mask=mask,
            cbar=True,
            linewidths=0.1,
            linecolor="white",
        )  # ,annot=annot,fmt=""#,mask=mask,

        g2.set_yticklabels(g2.get_yticklabels(), rotation=0, fontsize=6)
        g2.set_xticklabels(g2.get_xticklabels(), rotation=90, fontsize=6)
        axes.set_xlabel("")
        axes.set_ylabel("")

        fig.savefig(output_folder + stat_name + "_comparison_heatmap.pdf", format="pdf")
        # fig.savefig(output_folder+"fig_comparison_heatmap.jpeg", format='jpeg')
    else:


        fig, axes = plt.subplots(1, 1, figsize=(10, 10))
        fig.subplots_adjust(left=0.3, right=0.99, bottom=0.3, top=0.99)

        pivot_table = pivot_table.fillna(0)
        mask = pivot_table.T == 0
        print(mask.shape)
        print(pivot_table.T)

        g2 = sns.heatmap(
            pivot_table.T,
            cmap=palette,
            ax=axes,
            square=True,
            xticklabels=1,
            yticklabels=1,
            mask=mask,
            cbar=True,
            linewidths=0.1,
            linecolor="white",
        )  # ,annot=annot,fmt=""#,mask=mask,

        g2.set_yticklabels(g2.get_yticklabels(), rotation=0, fontsize=6)
        g2.set_xticklabels(g2.get_xticklabels(), rotation=90, fontsize=6)
        axes.set_xlabel("")
        axes.set_ylabel("")

        fig.savefig(output_folder + stat_name + "_comparison_heatmap.pdf", format="pdf")
        # fig.savefig(output_folder+"fig_comparison_heatmap.jpeg", format='jpeg')


def paint_final_table(final_table, output_file):

    out.apply_multiple_testing_correction(
        final_table, pval_col="empirical_pvalue", method="fdr_bh", threshold=0.05
    )

    with open(final_table) as f:
        table = pd.read_csv(f)

    table["log"] = -np.log10(table["empirical_pvalue"].values)
    table["log"] = table["log"].replace(
        to_replace=np.inf, value=np.max(table[table["log"] != np.inf]["log"]) + 0.5
    )

    setnames = list(set(table.setname.values.tolist()))

    f, ax = plt.subplots(1, figsize=(8, 20))
    g = sns.barplot(
        x="log",
        y="setname",
        hue="analysis",
        data=table,
        palette=sns.color_palette("muted", len(setnames)),
    )
    ax.axvline(
        x=1.3,
        ymin=0,
        ymax=len(setnames),
        alpha=0.5,
        color="k",
        linestyle="--",
        linewidth=0.5,
    )
    ax.axvline(
        x=4,
        ymin=0,
        ymax=len(setnames),
        alpha=0.5,
        color="k",
        linestyle="--",
        linewidth=0.5,
    )
    plt.subplots_adjust(left=0.4)
    plt.subplots_adjust(right=0.99)
    ax.set_xlabel("-log(empirical pvalue)")
    plt.savefig(output_file + ".pdf", f="pdf")
    plt.savefig(output_file + ".png", f="png")

    diz = {}
    analysis = []
    for a, tab in table.groupby(["analysis"]):
        analysis.append(a)
        tab = tab.sort_values(by=["empirical_pvalue", "observed"], ascending=False)
        tab = tab.reset_index()
        tab["rank"] = tab.index.values
        print(tab.head())
        diz[a] = {}
        tab = tab.sort_values(by=["setname"], ascending=False)
        print(tab.head())
        diz[a]["ranking"] = [
            (i) for i in zip(tab.index.values.tolist(), tab["setname"].values.tolist())
        ]
        diz[a]["tab"] = tab
        print(tab[["setname", "rank"]])

    M_corr = np.zeros((len(analysis), len(analysis)))
    M_pvalue = np.zeros((len(analysis), len(analysis)))
    for i in range(len(analysis)):
        for j in range(len(analysis)):
            s_corr, s_pvalue = stats.spearmanr(
                diz[analysis[i]]["tab"]["rank"].values,
                diz[analysis[j]]["tab"]["rank"].values,
                axis=0,
                nan_policy="propagate",
            )
            print(
                "%s, %s : spearman r--- corr: %f, p %f"
                % (analysis[i], analysis[j], s_corr, s_pvalue)
            )
            M_corr[i, j] = s_corr
            M_pvalue[i, j] = s_pvalue

    labels = [""]
    for k in analysis:
        labels.append(k)

    f, ax = plt.subplots(1, 2, figsize=(10, 5))
    plt.xticks(np.arange(6))
    plt.yticks(np.arange(6))
    ax[0].imshow(M_corr)
    ax[0].set_yticklabels(labels)
    ax[0].set_xticklabels(labels, rotation=90)
    ax[0].set_title("spearman r")

    ax[1].imshow(M_pvalue)
    ax[1].set_xticklabels(analysis, rotation=90)
    ax[1].set_title("pvalue")
    plt.savefig(output_file + "_corr.pdf", f="pdf")


def plot_adjacency(
    network: "network_filename",
    output_folder: "output_folder",
    prefix: "prefix for the file",
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

    print(nodelabels.shape)

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
        # plt.imshow(nx.adjacency_matrix(graph, nodelist=nodelist).toarray(), cmap='Greys')
        # plt.title("Adjacency Matrix")
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

    plt.savefig(output_folder + prefix + "_adjacency_matrix.png", f="png")
    plt.savefig(output_folder + prefix + "_adjacency_matrix.pdf", f="pdf")

    # fig,axes=plt.subplots(1, 3, figsize=(8,10),gridspec_kw={'width_ratios': [1,3,4],"wspace":0.025})
    # fig.subplots_adjust(left=0.2,right=0.80)

