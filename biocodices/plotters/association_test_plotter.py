import numpy as np
import seaborn as sns

from biocodices.plotters.base_plotter import BasePlotter
from biocodices.helpers.config import Config


class AssociationTestPlotter(BasePlotter):
    def __init__(self, association_test, plots_dir=None):
        self.assoc = association_test
        self.base_dir = plots_dir
        if plots_dir is None:
            self.base_dir = self.assoc.dataset.dir

    def draw_ax(self, ax, tests_to_plot, significance=0.05):
        df = self.assoc.result.sort_values(by=['CHR', 'BP'])
        palette_name = Config('plots')['assoc']['palette_name']
        colors = sns.color_palette(palette_name, len(tests_to_plot))
        for test in tests_to_plot:
            marker = '.' if test == 'UNADJ' else 'o'
            df[test].plot(ax=ax, marker=marker, linestyle='', logy=True,
                          label=test, color=colors.pop(0), linewidth=1)
        ax.invert_yaxis()
        ax.axhline(significance, linestyle='solid', color='green',
                   linewidth=1)
        ax.text(1, significance, '$P = {}$'.format(significance),
                color='green')
        ax.set_xlabel('')
        ax.set_xticks(np.arange(len(df.index)))
        ax.set_xticklabels(df.index, rotation=90)
        sns.despine(left=True)
        ax.legend()
        return ax
