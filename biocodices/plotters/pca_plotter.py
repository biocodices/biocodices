import seaborn as sns

from biocodices.plotters.base_plotter import BasePlotter
from biocodices.helpers.helpers import hide_spines_and_ticks, marker_from_sexcode
from biocodices.helpers.config import Config


class PCAPlotter(BasePlotter):
    def __init__(self, pca, plots_dir=None):
        """
        Expects a PCA object with 'results' and 'explained_variance'
        """
        self.pca = pca
        self.colors = Config('colors')  # FIXME use super()__init__()!
        self.base_dir = plots_dir  # FIXME should be in super too
        if plots_dir is None:
            self.base_dir = self.pca.dataset.dir
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

        by_phenotype = selected_components.groupby(level='phenotype')
        for phenotype, df in by_phenotype:
            color = self._new_color()
            for sexcode, values in df.groupby(level='sexcode'):
                marker = marker_from_sexcode(sexcode)
                kwargs = self._plot_kwargs(phenotype)
                x, y = components_to_plot
                label = 'Phenotype {}'.format(phenotype)
                values.plot(kind='scatter', x=x, y=y, ax=ax, label=label,
                            marker=marker, color=color, **kwargs)

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

    def _plot_kwargs(self, phenotype):
        kwargs = {
            # Generate a new color for a phenotype if there's no color defined
            # in the settings yml.
            # 'color': self.colors.get(phenotype, self._new_color()),
            # 'marker': self.plot_settings[phenotype]['marker'],
            'lw': self.plot_settings[phenotype]['linewidth'],
            'alpha': self.plot_settings[phenotype]['alpha'],
            's': self.plot_settings[phenotype]['markersize'],
            'zorder': self.plot_settings[phenotype]['zorder'],
        }
        return kwargs

    def _new_color(self):
        if not hasattr(self, '_more_colors'):
            palette_name = self.colors['QualitativePalette']
            phenotypes = self.pca.result.index.get_level_values('phenotype')
            number_of_phenotypes = len(phenotypes.unique())
            self._more_colors = sns.color_palette(palette_name,
                                                  number_of_phenotypes)
        return self._more_colors.pop(0)
