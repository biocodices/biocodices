import seaborn as sns
import ternary
from helpers.config import Config
from helpers.helpers import hide_spines_and_ticks


class AdmixturePlotter:
    def __init__(self, admixture, plots_dir=None):
        """
        Expects an admixture object that responds to #result
        """
        self.admixture = admixture
        self.colors = Config('colors')  # FIXME use super()__init__()!
        self.base_dir = plots_dir  # FIXME should be in super too
        if plots_dir is None:
            self.base_dir = self.admixture.dataset.source.plots_dir
        self.plot_settings = Config('plots')['admixture']

    def draw_ax(self, ax, K_to_plot, population_means=False):
        ancestries = self.admixture.result[K_to_plot]
        ancestries = self._reorder_samples_and_parse(ancestries)
        if population_means:
            ancestries = ancestries.groupby(level='population').mean()
            population_order = self.plot_settings['population_order']
            ancestries = ancestries.loc[population_order].dropna()

        plot_title = self._make_title(K_to_plot)
        palette = self._generate_palette(ancestries.columns)

        width = 0.5 if population_means else 1
        ancestries.plot(ax=ax, kind="bar", stacked=True, linewidth=0,
                        width=width, color=palette)

        self._plot_aesthetics(ax, plot_title, population_means, ancestries)

        if population_means:
            sns.despine(top=True, left=True, right=True)
        return ax

    def draw_triangle_ax(self, ax):
        ancestries = self.admixture.result[3]
        fig, tax = ternary.figure(scale=1, ax=ax)

        ancestries = ancestries[["EUR", "AFR", "AMR"]].dropna()
        by_population = ancestries.groupby(level="population", sort=False)
        for population, df in by_population:
            tax.scatter(df.values, label=population, s=45, marker='o',
                        color=self.colors[population])

        plot_title = self._make_title(3).replace(' - ', '\n')
        self._ternary_plot_aesthetics(tax, plot_title, ancestries)
        return tax

    def draw_cv_error(self, ax):
        Ks = list(self.admixture.cv_error.keys())
        ymin = self.admixture.cv_error.min()
        xmin = self.admixture.cv_error.idxmin()
        self.admixture.cv_error.plot(ax=ax, zorder=0, lw=0.75, color='grey',
                                     marker='o')
        ax.scatter(xmin, ymin, marker='o', zorder=1, s=85, edgecolor='black',
                   color='tomato', lw=1)
        ax.axvline(xmin, color='grey', linestyle='dotted')
        sns.despine(left=True, bottom=True)
        ax.set_ylabel("CV Error")
        ax.set_xlim(Ks[0] - 0.10, Ks[-1] + 0.10)
        ax.set_xticks(Ks)
        return ax

    def _reorder_samples_and_parse(self, ancestries):
        ancestries = ancestries.reset_index()
        # Use the same sample order through plots once defined.
        if not hasattr(self, 'plot_samples_order'):
            if 'AMR' in ancestries.columns:
                ancestries = ancestries.sort_values(['population', 'AMR'],
                                                    ascending=False)
            self.plot_samples_order = ancestries['sample']
        ancestries = ancestries.set_index('sample')
        ancestries = ancestries.loc[self.plot_samples_order]
        ancestries = ancestries.reset_index()
        ancestries = ancestries.drop(['region', 'family', 'sample'], axis=1)
        ancestries = ancestries.set_index("population")
        population_order = self.plot_settings['population_order']
        ancestries = ancestries.loc[population_order].dropna()
        return ancestries

    def _generate_palette(self, ancestry_labels):
        defined_colors = Config('colors')
        palette = []
        for ancestry in ancestry_labels:
            if ancestry in defined_colors:
                palette.append(defined_colors[ancestry])
        remaining = len(ancestry_labels) - len(palette)
        quali_palette = Config('colors')['QualitativePalette']
        palette += sns.color_palette(quali_palette, remaining)
        return palette

    def _plot_aesthetics(self, ax, plot_title, population_means, ancestries):
        ax.set_title(plot_title, y=1.01, family='serif')
        ax.legend_.set_visible(False)

        if not population_means:
            # Place the population labels in the middle of its samples
            population_order = ancestries.index.unique()
            N_by_population = ancestries.index.value_counts()[population_order]
            xlabels = N_by_population.cumsum() - N_by_population / 2
            ax.set_xticklabels(xlabels.index)
            ax.set_xticks(xlabels.values)

        ax.set_xlabel("")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
        ax.set_ylabel("Ancestría")
        ax.set_ylim([0, 1])
        ax.set_yticks([0, 1])
        #  ax.set_yticklabels([0, 1])
        #  hide_spines_and_ticks(ax, spines='all')

    def _ternary_plot_aesthetics(self, tax, title, ancestries):
        hide_spines_and_ticks(tax.get_axes(), spines="all")
        tax.boundary(linewidth=0.25)
        tax.clear_matplotlib_ticks()
        tax.set_title(title, position=(0.5, 1.15))
        tax.legend(frameon=False,  scatterpoints=1, bbox_to_anchor=(0.975, 1.125))
        tax.bottom_axis_label(ancestries.columns[0], position=(1, 0, 0), rotation=0)
        tax.right_axis_label(ancestries.columns[1], position=(-0.1, 1.2, -0.1), rotation=0)
        tax.left_axis_label(ancestries.columns[2], position=(0, 0, 1), rotation=0)

        return tax

    def _new_color(self):
        if not hasattr(self, '_more_colors'):
            palette_name = self.colors['QualitativePalette']
            regions = self.admixture.result.index.get_level_values('region')
            number_of_regions = len(regions.unique())
            self._more_colors = sns.color_palette(palette_name,
                                                  number_of_regions)
        return self._more_colors.pop(0)

    def _make_title(self, K):
        dataset = self.admixture.dataset
        return "{} - {} (K = {})".format(dataset.samplegroup.name,
                                         dataset.panel.name, K)
