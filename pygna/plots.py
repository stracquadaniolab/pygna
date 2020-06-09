from matplotlib.offsetbox import AnchoredText
from abc import ABC
import textwrap
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from palettable.colorbrewer.diverging import *
from palettable.colorbrewer.sequential import *


class PygnaFigure(ABC):
    """
    Abstract class that implements a general figure in Pygna. It has a class attribute representing the default ratios
    for the figures
    """
    ratios = {
        "volcano": [3, 4],
        "heatmap": [4, 3],
        "distribution": [4, 3],
        "barplot": [3, 4]
    }

    def __init__(self):
        super(PygnaFigure, self).__init__()
        self.filename = None

    def _save_fig(self) -> None:
        """
        This method saves the figure using the matplotlib framework
        """
        plt.savefig(self.filename)

    def _get_dimensions(self, fig_type: str, size: int) -> [int, int, int, int]:
        """
        This methods maps the ratios saved in the class with the figure type and returns the correct ration for each
        figure

        :param fig_type:  the figure type to be printed
        :param size:  the size of the plot
        :return: the width, height, fontsize, scalar of the figure
        """
        width = self.ratios[fig_type][0] * size
        height = self.ratios[fig_type][1] * size
        fontsize = size * 5
        if size <= 3:
            scalar = size ** 4
        else:
            scalar = size ** 3
        return int(width), int(height), int(fontsize), int(scalar)


class VolcanoPlot(PygnaFigure):
    """
    This class represent a Volcano Plot. It saves the value in the dataframe on a file.
    """

    def __init__(self, df: pd.DataFrame, output_file: str, loc: int = 2, p_col: str = "empirical_pvalue",
                 id_col: str = "setname_B", plotting_col: str = "observed", x_threshold: float = 0.1,
                 y_threshold: float = 0.1, y_label: str = "-log10(pvalue)", x_label: str = "z-score",
                 annotate: bool = False, size: int = 2):
        """
        :param df: the dataframe containing the data
        :param output_file: the output file
        :param loc: the location where to print the legend
        :param p_col: the column containing the empirical p-value values
        :param id_col: the column containing the setname
        :param plotting_col: the column to plot
        :param x_threshold: value of the threshold for the x axis
        :param y_threshold: value of the threshold for the y axis
        :param y_label: label for the y axis
        :param x_label: label for the x axis
        :param annotate: whether should be printed the annotaion table
        :param size: the size of the plots
        """
        super().__init__()
        self.df = df
        self.filename = output_file
        self.p_col = p_col
        self.id_col = id_col
        self.plotting_col = plotting_col
        self.x_thresh = x_threshold
        self.y_thresh = y_threshold
        self.y_label = y_label
        self.x_label = x_label
        self.annotate = annotate
        self.location = loc
        self.size = int(size)

        logging.info('Significant are considered when  %s > %f and %s > %f' % (self.p_col, self.y_thresh,
                                                                               self.plotting_col, self.x_thresh))
        self.df = self.df.sort_values(by=[p_col, plotting_col], ascending=True)
        not_sig = self._elaborate_not_sig_genes(self.df, self.plotting_col, self.p_col, self.y_thresh, self.x_thresh)
        sig = self._elaborate_sig_genes(self.df, self.plotting_col, self.p_col, self.y_thresh, self.x_thresh)
        text = [['Top 5 terms']]
        width, height, font_ratio, marker_scalar = super()._get_dimensions("volcano", self.size)
        fig, ax = plt.subplots(1, figsize=(width, height))

        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(font_ratio)
        ax.scatter(not_sig[plotting_col], not_sig[p_col], marker="o", alpha=1, edgecolors=None,
                   color=RdBu_4.mpl_colors[3], s=marker_scalar)
        if len(sig) < 1:
            logging.error("There are no significant terms, not plotting the volcano plot")
        else:
            ax.scatter(sig[plotting_col], sig[p_col], marker="+", alpha=1, edgecolors=None, color=RdBu_4.mpl_colors[0],
                       s=marker_scalar)
            if len(sig) > 5:
                ax.scatter(sig.iloc[:5, :][plotting_col], sig.iloc[:5, :][p_col], marker="*", alpha=1, s=marker_scalar,
                           edgecolors=None, color=RdBu_4.mpl_colors[0])
                self._append_name(sig.iloc[:5, :].iterrows(), self.id_col, text)

        ax.axhline(y=self.y_thresh, xmin=0, xmax=1, alpha=0.5, color='k', linestyle='--', linewidth=0.5)
        ax.axvline(x=self.x_thresh, ymin=0, ymax=1, alpha=0.5, color='k', linestyle='--', linewidth=0.5)
        ax.set_xlabel(self.x_label, fontsize=font_ratio)
        ax.set_ylabel(self.y_label, fontsize=font_ratio)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if self.annotate:
            texts = [sublist if type(sublist) == str else item for sublist in text for item in sublist]
            anchored_text = AnchoredText('\n'.join(texts), self.location, prop={'fontsize': font_ratio})
            ax.add_artist(anchored_text)

        super()._save_fig()

    def _elaborate_not_sig_genes(self, df, plotting_column, pvalue_col, y_thresh, x_thresh):
        """
        This method elaborates the genes which are not significant

        :param df: pd.datafrme, the dataframe to elaborate
        :param plotting_column: str, the column to apply the x threshold
        :param pvalue_col: str, the column to apply the y threshold
        :param y_thresh: float, the value of the x threshold
        :param x_thresh: float, the value of y threshold
        :return: pd.dataframe, the dataframe with the not significant genes
        """
        dataframe = df[(np.abs(df[plotting_column]) < x_thresh) | (df[pvalue_col] < y_thresh)].copy()
        return dataframe

    def _elaborate_sig_genes(self, df, plotting_column, pvalue_col, y_thresh, x_thresh):
        """
        This method elaborates the genes which are significant

        :param df: pd.datafrme, the dataframe to elaborate
        :param plotting_column: str, the column to apply the x threshold
        :param pvalue_col: str, the column to apply the y threshold
        :param y_thresh: float, the value of the x threshold
        :param x_thresh: float, the value of y threshold
        :return: pd.dataframe, the dataframe with the significant genes
        """
        dataframe = df[(np.abs(df[plotting_column]) >= x_thresh) & (df[pvalue_col] >= y_thresh)].copy()
        if len(dataframe) < 1:
            logging.warning("There are no significant terms")
        return dataframe

    def _append_name(self, dataset, id_col, text_list):
        """
        This method add the names to the legend of the plot

        :param dataset: df.dataframe, the dataset where are stored the names
        :param id_col: str, the column where the names are stored
        :param text_list: list, list with the names already stored
        :return: list, list with the names
        """
        for k, row in dataset:
            w = textwrap.wrap('-' + row[id_col].replace("_", " "), 30, break_long_words=False)
            text_list.append(w)
        return text_list
