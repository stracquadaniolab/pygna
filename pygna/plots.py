import matplotlib
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

    def _save_fig(self, filename):
        # TODO check the extension and save accordingly with the "format=" parameter
        plt.savefig(filename)

    def _get_dimensions(self, fig_type, size):
        width = self.ratios[fig_type][0]*size
        height = self.ratios[fig_type][1]*size
        fontsize = size*5
        scalar = size*3
        return width, height, fontsize, scalar


class VolcanoPlot(PygnaFigure):
    """
    This class represent a Volcano Plot
    """

    def __init__(self, df, output_file, p_col="empirical_pvalue", id_col="setname_B", plotting_col="observed",
                 x_threshold=0.1, y_threshold=0.1, y_label="-log10(pvalue)", x_label="z-score", annotate=False,
                 size=2):
        super().__init__()
        self.df = df
        self.output_file = output_file
        self.p_col = p_col
        self.id_col = id_col
        self.plotting_col = plotting_col
        self.x_thresh = x_threshold
        self.y_thresh = y_threshold
        self.y_label = y_label
        self.x_label = x_label
        self.annotate = annotate
        self.size = size

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

        super()._save_fig(self.output_file)

    def _elaborate_not_sig_genes(self, df, plotting_column, pvalue_col, y_thresh, x_thresh):
        dataframe = df[(df[plotting_column] < x_thresh) | (df[pvalue_col] < y_thresh)].copy()
        return dataframe

    def _elaborate_sig_genes(self, df, plotting_column, pvalue_col, y_thresh, x_thresh):
        dataframe = df[(np.abs(df[plotting_column]) >= x_thresh) & (df[pvalue_col] >= y_thresh)].copy()
        if len(dataframe) < 1:
            logging.info('there are no significant terms')
        return dataframe

    def _append_name(self, condition, id_col, text_list):
        for k, row in condition:
            w = textwrap.wrap('-' + row[id_col].replace("_", " "), 30, break_long_words=False)
            text_list.append(w)
        return text_list
