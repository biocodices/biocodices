from os.path import join
from math import ceil

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


class VariantCallingMetricsPlotter:
    def __init__(self, df):
        self.data = df

    def plot_and_savefig(self, out_dir):
        sns.set_context('notebook')
        sns.set_style('white')

        #  plot_w = 3 + len(self.data['sample'].unique())
        #  plot_h = 3.5
        #  plots_per_row = 3

        #  n_plots = len(self.data.columns) - 2
        #  n_rows = ceil(n_plots / plots_per_row)
        #  n_cols = ceil(n_plots / n_rows)
        #  ax_ids = list(np.arange(n_plots) + 1)

        #  fig = plt.figure()
        #  fig.set_figheight(plot_h * n_rows)
        #  fig.set_figwidth(plot_w * n_cols)

        #  for i, category in enumerate(self.data.columns):
            #  if category in ['CATEGORY', 'sample']:
                #  continue

            #  ax = fig.add_subplot(n_rows, n_cols, ax_ids.pop(0))
            #  self.draw_ax(ax, category)
            #  if i == 0:
                #  ax.legend()
            #  else:
                #  ax.legend_.set_visible(False)

        plt.tight_layout()
        fp = join(out_dir, 'alignment_metrics')
        plt.savefig(fp, dpi=300, bbox_inches='tight')

    def draw_ax(self, ax, category):
        #  sns.barplot(ax=ax, x='sample', y=category, hue='CATEGORY',
                    #  data=self.data)
        #  ax.set_xlabel('')
        #  ax.set_ylabel('')
        #  ax.set_title(category)
        #  if 'PCT' in category:
            #  ax.set_ylim([0, 1])
        sns.despine(left=True)
