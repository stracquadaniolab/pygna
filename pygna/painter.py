import logging
import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pygna.output as out
import pygna.parser as ps
import pygna.plots as pp
import seaborn as sns
from palettable.colorbrewer.diverging import *
from palettable.colorbrewer.sequential import *


# TODO Refactor and check this
def volcano_plot(df, output_file, p_col, id_col, plotting_col, threshold_x, threshold_y, y_label, x_label, annot, loc):
    pp.VolcanoPlot(df=df, output_file=output_file, p_col=p_col, id_col=id_col, plotting_col=plotting_col,
                   x_threshold=threshold_x, y_threshold=threshold_y, y_label=y_label, x_label=x_label, annotate=annot,
                   loc=loc)


#######################################################
######## COMMAND FUNCTIONS ############################
########################################################


def paint_datasets_stats(table_filename: 'pygna results table',
                         output_file: 'figure file, use pdf or png extension',
                         alternative="greater"
                         ):
    '''
    This function plots the results of of a GNT test.
    Pass the results table generated by one of the functions
    and the output figure file (png or pdf).
    In case you are using a SP test, pass also 'less' as
    an alternative.
    '''

    palette_binary = RdBu_4.mpl_colors[0::3]

    with open(table_filename, "r") as f:
        table = pd.read_csv(f, sep=",")

    stat_name = table["analysis"][0]
    n_permutations = table["number_of_permutations"][0]
    data = pd.DataFrame()

    # If greater the test is between 0 and + inf
    # else is bwteen -inf and 0, I need to set the right bounds
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
        plt.savefig(output_file + '.png', format="png")


def paint_comparison_matrix(table_filename: 'pygna comparison output',
                            output_file: 'output figure file, specify png or pdf file',
                            rwr: 'use rwr is the table comes from a rwr analysis' = False,
                            single_geneset: 'use true if the comparison has been done for a single file' = False,
                            annotate: 'set true if uou want to print the pvalue inside the cell' = False):
    '''
    This function plots the results of of a GNA test.
    Pass the results table generated by one of the functions and the output figure file (png or pdf).
    With rwr you can specify whether the test is a rwr association, in this case a different palette and limits are sets.
    Specify if the results are obtained using association with only one genesets (multiple setnames in the same file).
    Pass the annotate flag to have the pvalue annotation on the plot
    '''

    # If rwr we want sequential, non divergent
    if rwr:
        palette = OrRd_9.mpl_colors
    else:
        palette = RdBu_11.mpl_colors

    with open(table_filename, "r") as f:
        table = pd.read_csv(f, sep=",")

    stat_name = table["analysis"][0]
    n_permutations = table["number_of_permutations"][0]

    # If we only have used one geneset, we need to add 1 col and 1 row
    # To have a symmetrical matrix for the pivot
    if single_geneset:
        table = table.loc[:, ['observed', 'setname_A', 'setname_B', 'empirical_pvalue']]
        for i in set(table['setname_A'].values.tolist()).union(set(table['setname_B'].values.tolist())):
            table = table.append({'setname_A': i, 'setname_B': i, 'observed': 0, 'empirical_pvalue': 1},
                                 ignore_index=True)

    pivot_table = table.pivot(values="observed", index="setname_A", columns="setname_B")
    pivot_table = pivot_table.fillna(0)

    # we make sure we have a triangual matrix (if I have A/B and not B/A
    # the pivot is leaving some cells NA )
    matrix = (pivot_table.T + pivot_table) - pivot_table.T * np.eye(len(pivot_table))
    if single_geneset:
        mask = np.triu(np.ones((len(pivot_table), len(pivot_table))))
    else:
        mask = np.ones((len(pivot_table), len(pivot_table))) - np.tril(np.ones((len(pivot_table), len(pivot_table))))

    # Create annotation from pvalues
    if annotate:
        annot = table.pivot(
            values="empirical_pvalue", index="setname_A", columns="setname_B"
        )
        annot = annot.fillna(0)
        annot = (annot.T.values + annot.values) - annot.T.values * np.eye(len(annot))
    else:
        annot = False

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
        fig.savefig(output_file + '.png', format="png")


def paint_volcano_plot(table_filename: 'pygna comparison output',
                       output_file: 'output figure file, specify png or pdf file',
                       rwr: 'use rwr is the table comes from a rwr analysis' = False,
                       id_col="setname_B",
                       threshold_x=0,
                       threshold_y=2,
                       annotate=False):
    '''
    This function plots the results of of a GNA test of association
    of a single geneset against multiple pathways.
    Pass the results table generated by one of the functions
    and the output figure file (png or pdf).
    From the results table, a multiple testing correction is applied
    and the results are those plotted.
    The defined threshold are for x: zscore and y: -log10(pvalue)

    '''

    out.apply_multiple_testing_correction(
        table_filename, pval_col="empirical_pvalue", method="fdr_bh", threshold=0.1
    )

    with open(table_filename, "r") as f:
        df = pd.read_csv(f, sep=",")

    stat_name = df["analysis"][0]
    n_permutations = df["number_of_permutations"][0]

    # Normalise the plotting value
    df['zscore'] = (df['observed'] - df['mean(null)']) / np.sqrt(df['var(null)'].values)

    # If it's a shortest path, mirror the zscore
    loc = 2
    if not rwr:
        loc = 1

    # When pvalue==0 the -log would be infinite, hence we replace
    # pvalue with the permutation resolution - a tenth of the resolution
    # So that the plotted value is going to be 1 + max_log
    sig_th = 1 / n_permutations - 1 / n_permutations / 10
    df['bh_pvalue'] = df['bh_pvalue'] + (df['bh_pvalue'] == 0.0) * sig_th

    # transform in -log10(pvalue)
    df['-log10(p)'] = -np.log10(df['bh_pvalue'].values)

    volcano_plot(df, output_file, p_col='-log10(p)', id_col=id_col, plotting_col="zscore", threshold_x=threshold_x,
                 threshold_y=threshold_y, y_label='-log10(pvalue)', x_label='z-score', annot=annotate, loc=loc)


#######################################################
#################### EXTRA ############################
########################################################


def paint_volcano_plot(table_filename: 'pygna comparison output',
                       output_file: 'output figure file, specify png or pdf file',
                       rwr: 'use rwr is the table comes from a rwr analysis' = False,
                       id_col="setname_B",
                       threshold_x=0,
                       threshold_y=2,
                       annotate=False):
    '''
    This function plots the results of of a GNA test of association
    of a single geneset against multiple pathways.
    Pass the results table generated by one of the functions
    and the output figure file (png or pdf).

    From the results table, a multiple testing correction is applied
    and the results are those plotted.

    The defined threshold are for x: zscore and y: -log10(pvalue)

    '''

    out.apply_multiple_testing_correction(
        table_filename, pval_col="empirical_pvalue", method="fdr_bh", threshold=0.1
    )

    with open(table_filename, "r") as f:
        df = pd.read_csv(f, sep=",")

    stat_name = df["analysis"][0]
    n_permutations = df["number_of_permutations"][0]

    # Normalise the plotting value
    df['zscore'] = (df['observed'] - df['mean(null)']) / np.sqrt(df['var(null)'].values)

    # If it's a shortest path, mirror the zscore
    loc = 2
    if not rwr:
        loc = 1

    # When pvalue==0 the -log would be infinite, hence we replace
    # pvalue with the permutation resolution - a tenth of the resolution
    # So that the plotted value is going to be 1 + max_log
    sig_th = 1 / n_permutations - 1 / n_permutations / 10
    df['bh_pvalue'] = df['bh_pvalue'] + (df['bh_pvalue'] == 0.0) * sig_th

    # transform in -log10(pvalue)
    df['-log10(p)'] = -np.log10(df['bh_pvalue'].values)

    volcano_plot(df, output_file, p_col='-log10(p)', id_col=id_col, plotting_col="zscore", threshold_x=threshold_x,
                 threshold_y=threshold_y, y_label='-log10(pvalue)', x_label='z-score', annot=annotate, loc=loc)


#######################################################
#################### EXTRA ############################
########################################################


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
    This function has been mostly used for the generation of plots in the paper.
    Please raise an issue if you want it to be improved.
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
        matrix = matrix[0: int(size), 0: int(size)]
        nodelabels = nodelabels[:, 0: int(size)]

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
        f.savefig(output_file + '.png', format="png")


def stars(pvalue) -> str:
    s = ""
    if pvalue > 0.05:
        s = "ns"
    elif 0.01 < pvalue <= 0.05:
        s = "*"
    elif 0.001 < pvalue <= 0.01:
        s = "**"
    elif 0.0001 < pvalue <= 0.001:
        s = "***"
    else:
        s = "****"
    return s
