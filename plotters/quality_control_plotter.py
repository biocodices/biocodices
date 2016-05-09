#  import numpy as np
#  import seaborn as sns

from plotters.base_plotter import BasePlotter
#  from helpers.config import Config


class QualityControlPlotter(BasePlotter):
    def __init__(self, quality_control, plots_dir=None):
        self.quali = quality_control
        self.base_dir = plots_dir or self.quali.dataset.dir

    def draw_ax(self, ax, tests_to_plot, significance=0.05):
        raise Exception('Not implemented yet')
