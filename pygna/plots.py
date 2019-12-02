from matplotlib.offsetbox import AnchoredText
from abc import ABC
import textwrap
import logging
import numpy as np
import matplotlib.pyplot as plt


class PygnaFigure(ABC):
    ratios = {
        "volcano": [3, 4],
        "heatmap": [4, 3],
        "distribution": [4, 3],
        "barplot": [3, 4]
    }

    def __init__(self):
        super(PygnaFigure, self).__init__()
        self.filename = None

    def _save_fig(self):
        """
        This method saves the figure using the matplotlib framework
        :return: null
        """
        # TODO check the extension and save accordingly with the "format=" parameter
        plt.savefig(self.filename)

    def _get_dimensions(self, fig_type, size):
        """
        This methods maps the ratios saved in the class with the figure type and returns the correct ration for each
        figure
        :param fig_type: str, the figure type to be printed
        :param size: int, the size of the plot
        :return: int, int, int, int the width, height, fontsize, scalar of the figure
        """
        width = self.ratios[fig_type][0]*size
        height = self.ratios[fig_type][1]*size
        fontsize = size*5
        if size <= 3:
            scalar = size**4
        else:
            scalar = size**3
        return int(width), int(height), int(fontsize), int(scalar)


class VolcanoPlot(PygnaFigure):
    """
    This class represent a Volcano Plot
    """

    def __init__(self, df, output_file, p_col="empirical_pvalue", id_col="setname_B", plotting_col="observed",
                 x_threshold=0.1, y_threshold=0.1, y_label="-log10(pvalue)", x_label="z-score", annotate=False,
                 size=2):
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
        self.size = int(size)

        logging.info('Significant are considered when  %s > %f and %s > %f' % (self.p_col, self.y_thresh,
                                                                               self.plotting_col, self.x_thresh))
        self.df = self.df.sort_values(by=[p_col, plotting_col], ascending=False)
        not_sig = self._elaborate_not_sig_genes(self.df, self.plotting_col, self.p_col, self.y_thresh, self.x_thresh)
        sig = self._elaborate_sig_genes(self.df, self.plotting_col, self.p_col, self.y_thresh, self.x_thresh)
        text = [['Top 5 terms']]
        width, height, font_ratio, marker_scalar = super()._get_dimensions("volcano", self.size)

        fig, ax = plt.subplots(1, figsize=(width, height))
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(font_ratio)
        ax.scatter(not_sig[plotting_col], not_sig[p_col], marker="o", alpha=1, edgecolors=None, color='blue',
                   s=marker_scalar)
        if len(sig) < 1:
            logging.error("There are no significant terms, not plotting the volcano plot")
        else:
            ax.scatter(sig[plotting_col], sig[p_col], marker="+", alpha=1, edgecolors=None, color='red',
                       s=marker_scalar)
            if len(sig) > 5:
                ax.scatter(sig.iloc[:5, :][plotting_col], sig.iloc[:5, :][p_col], marker="*", alpha=1, s=marker_scalar,
                           edgecolors=None, color='red')
                self._append_name(sig.iloc[:5, :].iterrows(), self.id_col, text)

        ax.axhline(y=self.y_thresh, xmin=0, xmax=1, alpha=0.5, color='k', linestyle='--', linewidth=0.5)
        ax.axvline(x=self.x_thresh, ymin=0, ymax=1, alpha=0.5, color='k', linestyle='--', linewidth=0.5)
        ax.set_xlabel(self.x_label, fontsize=font_ratio)
        ax.set_ylabel(self.y_label, fontsize=font_ratio)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if self.annotate:
            texts = [sublist if type(sublist) == str else item for sublist in text for item in sublist]
            anchored_text = AnchoredText('\n'.join(texts), loc=2, prop={'fontsize': font_ratio})
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
        dataframe = df[(df[plotting_column] < x_thresh) | (df[pvalue_col] < y_thresh)].copy()
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
