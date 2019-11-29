from abc import ABC, abstractmethod

import logging
import numpy as np
import matplotlib.pyplot as plt


class PygnaFigure(ABC):
    def __init__(self):
        super(PygnaFigure, self).__init__()


class VolcanoPlot(PygnaFigure):
    """
    This class represent a Volcano Plot
    """

    def __init__(self, df, output_file, p_col="empirical_pvalue", id_col="setname_B", plotting_col="observed",
                 x_threshold=0.1, y_threshold=0.1, y_label='-log10(pvalue)', x_label='z-score', annotate=False):
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
        nsig = self._elaborate_not_sig_genes(self.df, self.plotting_col, self.p_col, self.y_thresh, self.x_thresh)
        sig = self._elaborate_sig_genes(self.df, self.plotting_col, self.p_col, self.y_thresh, self.x_thresh)
        fig, ax = plt.subplots(1, figsize=(8, 10))
        ax.scatter(nsig[plotting_col], nsig[p_col], marker="o", s=20, alpha=1, edgecolors=None, color='blue')
        texts = []
        texts.append(['Top 5 terms'])




    def _elaborate_not_sig_genes(self, df, plotting_column, pvalue_col, y_thresh, x_thresh):
        dataframe = df[(df[plotting_column] < x_thresh) | (df[pvalue_col] < y_thresh)].copy()
        return dataframe

    def _elaborate_sig_genes(self, df, plotting_column, pvalue_col, y_thresh, x_thresh):
        dataframe = df[(np.abs(df[plotting_column]) >= x_thresh) & (df[pvalue_col] >= y_thresh)].copy()
        if len(dataframe) < 1:
            logging.info('there are no significant terms')
        return dataframe
