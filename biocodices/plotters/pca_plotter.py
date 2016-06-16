import seaborn as sns

from biocodices.plotters.base_plotter import BasePlotter
from biocodices.helpers.plotting import (hide_spines_and_ticks,
                                         marker_from_sexcode)
from biocodices.programs import Plink


class PCAPlotter(BasePlotter):
    def __init__(self, pca, plots_dir=None, palette=None):
        """
        Expects a PCA object with 'results' and 'explained_variance'
        """
        self.pca = pca
        self.base_dir = plots_dir  # FIXME should be in super too
        if plots_dir is None:
            self.base_dir = self.pca.dataset.dir
        self.explained_variance = self.pca.explained_variance
        self.palette = palette or 'husl'

    def draw_ax(self, ax, components_to_plot, show_ticks, **kwargs):
        """
        Draws a scatterplot of the first two columns in eigenvalues.
        Passees **kwargs to pandas .plot().
        """
        if len(components_to_plot) != 2:
            error_msg = 'I only know how to plot exactly TWO components. '
            error_msg += 'I received: {}'.format(components_to_plot)
            raise ValueError(error_msg)
        selected_components = self.pca.result[components_to_plot]

        by_phenotype = selected_components.groupby(level='phenotype')
        for phenotype, df in by_phenotype:
            color = self._new_color()
            print(phenotype, color)
            for sex_code, values in df.groupby(level='sexcode'):
                marker = marker_from_sexcode(sex_code)
                x, y = components_to_plot
                label = '{}, {}'.format(
                    Plink.phenotype_codes.get(phenotype, 'Case/Control Unknown'),
                    Plink.sex_codes.get(sex_code, 'Sex Unknown'),
                )
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

    def _new_color(self):
        if not hasattr(self, '_colors'):
            phenotypes = self.pca.result.index.get_level_values('phenotype')
            number_of_phenotypes = len(phenotypes.unique())
            self._colors = sns.color_palette(self.palette,
                                             number_of_phenotypes)
        return self._colors.pop(0)
