import matplotlib.pyplot as plt
import numpy as np
from math import sin, cos, degrees
from plotters.pca_plotter import PCAPlotter
from helpers.config import Config


class BasePCA:
    """
    Don't instantiate this class; use instead SmartPCA or SklearnPCA.
    """
    def plot(self, ax=None, components_to_plot=['PC1', 'PC2'], rotate=False,
             show_ticks=False):
        self.plotter = PCAPlotter(self, self.dataset.source.plots_dir)
        self.inverted_x = self.inverted_y = False
        self.invert_axes(components_to_plot)
        self.rotated = rotate
        if rotate:
            self.rotate_the_results(components_to_plot)

        if ax is None:
            _, ax = plt.subplots(figsize=(5, 5))
        self.plotter.draw_ax(ax, components_to_plot, show_ticks=show_ticks)
        return ax

    def __repr__(self):
        s = '<{} for {}>'.format(self.__class__.__name__,
                                 self.dataset.full_label)
        return s

    def savefig(self, filename=None):
        if filename is None:
            filename = '{}.{}'.format(self.dataset.label, type(self).__name__)
        self.plotter.savefig(filename)

    def write_result_csvs(self):
        self.result.to_csv(self._output_filepath('eigenvecs') + '.csv')
        self.explained_variance.to_csv(self._output_filepath('eigenvals') + '.csv')

    def _output_filepath(self, ext):
        return self.dataset.bedfile + '.' + ext

    def invert_axes(self, components):
        if len(components) != 2:
            raise ValueError('I need exactly two components as x and y.')

        x, y = self.result[components[0]], self.result[components[1]]
        x_domain, y_domain = [x.min(), x.max()], [y.min(), y.max()]
        x_mean, y_mean = np.mean(x_domain), np.mean(y_domain)

        ref_population = Config('plots')['PCA']['reference_top']
        x_refpop = x.xs(ref_population, level='population')
        y_refpop = y.xs(ref_population, level='population')
        x_median, y_median = np.median(x_refpop), np.median(y_refpop)

        ref_population_on_top = y_median > y_mean
        ref_population_in_the_left = x_median < x_mean

        if not ref_population_in_the_left:
            self.result[components[0]] = x.apply(lambda n: n * -1)
            self.inverted_x = True
        if not ref_population_on_top:
            self.result[components[1]] = y.apply(lambda n: n * -1)
            self.inverted_y = True

    def _define_angle_of_rotation(self, components):
        if len(components) != 2:
            raise ValueError('I need exactly two components as x and y.')

        reference_top = Config('plots')['PCA']['reference_top']
        reference_bottom = Config('plots')['PCA']['reference_bottom']
        top_samples = self.result.xs(reference_top, level='population')
        bottom_samples = self.result.xs(reference_bottom, level='population')
        top_leftmost = top_samples[components[0]].idxmin()
        bottom_leftmost = bottom_samples[components[0]].idxmin()
        #  top_x, top_y = top_samples[components].quantile(0.20)
        #  bottom_x, bottom_y = bottom_samples[components].quantile(0.20)
        top_x, top_y = top_samples.loc[top_leftmost][components]
        bottom_x, bottom_y = bottom_samples.loc[bottom_leftmost][components]
        vector = ((top_x - bottom_x), (top_y - bottom_y))
        vectors = [(0, 1), vector]
        unit_vectors = [v / np.linalg.norm(v) for v in vectors]
        angle = np.arccos(np.clip(np.dot(*unit_vectors), -1.0, 1.0))
        if top_x < bottom_x:
            angle *= -1
        return angle

    def rotate_the_results(self, components):
        if len(components) != 2:
            raise ValueError('I need exactly two components as x and y.')

        self.invert_axes(components)
        angle = self._define_angle_of_rotation(components)
        original_values = self.result[components]
        x_name, y_name = components
        rotated_values = original_values.copy()
        for index, x, y in original_values.itertuples():
            new_x = x * cos(angle) - y * sin(angle)
            new_y = y * cos(angle) + x * sin(angle)
            rotated_values.loc[index][x_name] = new_x
            rotated_values.loc[index][y_name] = new_y

        for component in rotated_values.columns:
            self.result[component] = rotated_values[component]

        self.rotated = True
        self.rotation_angle = degrees(angle)
