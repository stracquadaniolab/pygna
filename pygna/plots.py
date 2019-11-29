import textwrap
from abc import ABC, abstractmethod

import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText


class PygnaFigure(ABC):
    def __init__(self):
        super(PygnaFigure, self).__init__()

    def _save_fig(self, filename):
        # TODO check the extension and save accordingly with the "format=" parameter
        plt.savefig(filename)


class VolcanoPlot(PygnaFigure):
    """
    This class represent a Volcano Plot
    """

    def __init__(self, df, output_file, p_col="empirical_pvalue", id_col="setname_B", plotting_col="observed",
                 x_threshold=0.1, y_threshold=0.1, y_label="-log10(pvalue)", x_label="z-score", annotate=False):
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

        logging.info('Significant are considered when  %s > %f and %s > %f' % (self.p_col, self.y_thresh,
                                                                               self.plotting_col, self.x_thresh))
        not_sig = self._elaborate_not_sig_genes(self.df, self.plotting_col, self.p_col, self.y_thresh, self.x_thresh)
        sig = self._elaborate_sig_genes(self.df, self.plotting_col, self.p_col, self.y_thresh, self.x_thresh)
        fig, ax = plt.subplots(1, figsize=(8, 10))
        ax.scatter(not_sig[plotting_col], not_sig[p_col], marker="o", s=20, alpha=1, edgecolors=None, color='blue')
        text = [['Top 5 terms']]

        if len(sig) < 1:
            logging.error("There are no significant terms, not plotting the volcano plot")
        elif len(sig) > 5:
            ax.scatter(sig.iloc[:5, :][plotting_col], sig.iloc[:5, :][p_col], marker="*", s=70, alpha=1,
                       edgecolors=None, color='red')
            ax.scatter(sig.iloc[5:, :][plotting_col], sig.iloc[5:, :][p_col], marker="+", s=50, alpha=.5,
                       edgecolors=None, color='red')
            self._append_name(sig.iloc[:5, :].iterrows(), self.id_col, text)
        else:
            ax.scatter(sig[plotting_col], sig[p_col], marker="*", s=50, alpha=1, edgecolors=None, color='red')
            self._append_name(sig.iterrows(), self.id_col, text)

        ax.axhline(y=self.y_thresh, xmin=0, xmax=1, alpha=0.5, color='k', linestyle='--', linewidth=0.5)
        ax.axvline(x=self.x_thresh, ymin=0, ymax=1, alpha=0.5, color='k', linestyle='--', linewidth=0.5)
        ax.set_xlabel(self.x_label, fontsize=12)
        ax.set_ylabel(self.y_label, fontsize=12)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if self.annotate:
            texts = [sublist if type(sublist) == str else item for sublist in text for item in sublist]
            anchored_text = AnchoredText('\n'.join(texts), loc=2, prop={'fontsize': 12})
            ax.add_artist(anchored_text)

        PygnaFigure._save_fig(self, self.output_file)

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
