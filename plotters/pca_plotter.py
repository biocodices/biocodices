import seaborn as sns

from plotters.base_plotter import BasePlotter
from helpers.helpers import hide_spines_and_ticks
from helpers.config import Config


class PCAPlotter(BasePlotter):
    def __init__(self, pca, plots_dir=None):
        """
        Expects a PCA object with 'results' and 'explained_variance'
        """
        self.pca = pca
        self.colors = Config('colors')  # FIXME use super()__init__()!
        self.base_dir = plots_dir  # FIXME should be in super too
        if plots_dir is None:
            self.base_dir = self.pca.dataset.source.plots_dir
        self.plot_settings = Config('plots')['PCA']
        self.explained_variance = self.pca.explained_variance

    def draw_ax(self, ax, components_to_plot, show_ticks):
        """
        Draws a scatterplot of the first two columns in eigenvalues
        """
        if len(components_to_plot) != 2:
            error_msg = 'I only know how to plot exactly TWO components. '
            error_msg += 'I received: {}'.format(components_to_plot)
            raise ValueError(error_msg)
        selected_components = self.pca.result[components_to_plot]

        grouped_results = selected_components.groupby(level='population')
        for population, values in grouped_results:
            kwargs = self._plot_kwargs(population)
            x, y = components_to_plot
            values.plot(kind='scatter', x=x, y=y, ax=ax, label=population,
                        **kwargs)

        # Set the axes labels
        xlabel_prefix = '-' if self.pca.inverted_x else ''
        ylabel_prefix = '-' if self.pca.inverted_y else ''
        xlabel_suffix = ''
        if self.pca.rotated:
            xlabel_suffix = '\nRotated {}°'.format(int(self.pca.rotation_angle))

        xvariance = self.explained_variance.ix[x]['percentage']
        xlabel = "{}{}: {}{}".format(xlabel_prefix, x, xvariance,
                                     xlabel_suffix)
        ax.set_xlabel(xlabel)
        yvariance = self.explained_variance.ix[y]['percentage']
        ylabel = "{}{}: {}".format(ylabel_prefix, y, yvariance)
        ax.set_ylabel(ylabel)

        if not show_ticks:
            #  Remove non-data ink
            ax.tick_params(axis="x", bottom="off", top="off")
            ax.tick_params(axis="y", left="off", right="off")
            hide_spines_and_ticks(ax, 'all')

        return ax

    def _plot_kwargs(self, population):
        primary = population in self.plot_settings['primary_populations']
        importance = 'primary' if primary else 'secondary'
        kwargs = {
            # Generate a new color for a population if there's no color defined
            # in the settings yml.
            'color': self.colors.get(population, self._new_color()),
            'marker': self.plot_settings[importance]['marker'],
            'lw': self.plot_settings[importance]['linewidth'],
            'alpha': self.plot_settings[importance]['alpha'],
            's': self.plot_settings[importance]['markersize'],
            'zorder': self.plot_settings[importance]['zorder'],
        }
        return kwargs

    def _new_color(self):
        if not hasattr(self, '_more_colors'):
            palette_name = self.colors['QualitativePalette']
            populations = self.pca.result.index.get_level_values('population')
            number_of_populations = len(populations.unique())
            self._more_colors = sns.color_palette(palette_name,
                                                  number_of_populations)
        return self._more_colors.pop(0)
