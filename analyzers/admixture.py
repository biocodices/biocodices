import re
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

from os import chdir, getcwd
from os.path import expanduser, isfile
from plotters.admixture_plotter import AdmixturePlotter
from helpers.config import Config


class Admixture:
    _EXECUTABLE = expanduser(Config('executables')['admixture'])

    def __init__(self, dataset):
        self.dataset = dataset
        self.Pfiles = {}
        self.Qfiles = {}
        self.logfiles = {}
        self.result = {}
        self.cv_error = pd.Series([], index=[])
        self.cv_error.index.name = 'K'
        self.config = Config('plots')['admixture']
        self.plotter = AdmixturePlotter(self, self.dataset.source.plots_dir)

    def __repr__(self):
        return "<Admixture for {}, Ks={}>".format(self.dataset.full_label,
                                                  list(self.result.keys()))

    def run(self, Ks, cores, infer_components=False, overwrite=False):
        if not hasattr(Ks, '__iter__'):
            Ks = [Ks]

        for K in Ks:
            self.logfiles[K] = self._output_filepath(K, 'log')
            self.Pfiles[K] = self._output_filepath(K, 'P')
            self.Qfiles[K] = self._output_filepath(K, 'Q')
            analysis_exists = isfile(self.Pfiles[K]) or isfile(self.Qfiles[K])
            if not analysis_exists or overwrite:
                self._call_admixture(K, cores)
            self.result[K] = self._read_ancestry_file(K)
            self.cv_error.loc[K] = self._read_cv_error_file(K)

            regions = self.result[K].index.get_level_values('region').unique()
            if infer_components and len(regions) >= 3:
                self._assign_regions_to_clusters(self.result[K])
            self._infer_clusters_from_reference_population(self.result[K])
            self.result[K] = self._reorder_clusters(self.result[K])

    def population_means(self, K):
        if K not in self.result:
            raise ValueError("I don't have results for K={}.".format(K))
        return self.result[K].groupby(level='population').mean()

    def plot(self, K, population_means=False, ax=None):
        if ax is None:
            _, ax = plt.subplots(figsize=(15, 2.5))
        self.plotter.draw_ax(ax, K, population_means=population_means)
        return ax

    def plot_triangle(self, ax=None):
        if 3 not in self.result:
            raise Exception("I don't have results for K=3!")
        if ax is None:
            _, ax = plt.subplots(figsize=(7, 6))
        return self.plotter.draw_triangle_ax(ax=ax)

    def plot_cv_error(self, ax=None):
        if ax is None:
            _, ax = plt.subplots(figsize=(8, 4))
        return self.plotter.draw_cv_error(ax=ax)

    def savefig(self, filename=None):
        if filename is None:
            filename = '{}.{}'.format(self.dataset.label, type(self).__name__)
        self.plotter.savefig(filename)

    def _assign_regions_to_clusters(self, ancestries_df):
        means = ancestries_df.groupby(level='region').mean()
        max_region_per_component = means.idxmax(axis=0).to_dict()
        max_component_per_region = means.idxmax(axis=1).to_dict()

        # I accept an association of region-cluster if two conditions are met:
        # - The component must have the max avg. ancestry value for that
        #   region.
        # - For the region, that component must have the max avg. ancestry
        #   value.
        # I know it's kind of a tongue twister, but check the result DataFrames
        # to get what I'm saying. Max value column-wise must coincide with max
        # value row-wise, otherwise it would be a dubious guess.
        guesses = {}
        for region, guessed_component in max_component_per_region.items():
            if region == max_region_per_component[guessed_component]:
                guesses[guessed_component] = region

        ancestries_df.rename(columns=guesses, inplace=True)

    def _infer_clusters_from_reference_population(self, ancestries):
        reference_population = self.config['reference_population']
        reference_ancestries = self.config['reference_population_ancestries']
        by_population = ancestries.groupby(level='population')
        mean_ancestries = by_population.mean().loc[reference_population]
        components_order = mean_ancestries.sort_values(ascending=False).index

        guesses = {}
        for component, ancestry in zip(components_order, reference_ancestries):
            if type(component) != int:
                continue  # Don't re-guess already known ancestries
            if ancestry not in ancestries.columns:
                guesses[component] = ancestry

        if len(guesses) > 0:
            ancestries.rename(columns=guesses, inplace=True)

    def _reorder_clusters(self, ancestries):
        cluster_order = Config('plots')['admixture']['cluster_order']
        ordered = [c for c in cluster_order if c in ancestries.columns]
        the_rest = [c for c in ancestries.columns if c not in ordered]
        ancestries = ancestries[ordered + the_rest]
        return ancestries

    def _read_ancestry_file(self, K):
        result = pd.read_table(self.Qfiles[K], sep='\s+', header=None)
        result.columns = range(1, K+1)
        result.index = self.dataset.samplegroup.samples.index
        result = result.join(self.dataset.samplegroup.samples)
        multi_index = ['region', 'population', 'family', 'sample']
        result = result.reset_index().set_index(multi_index)
        return result

    def _read_cv_error_file(self, K):
        log_fn = self.logfiles[K]
        with open(log_fn, 'r') as logfile:
            cv_error_lines = [line.strip() for line in logfile.readlines()
                              if re.match('CV error', line)]

            if len(cv_error_lines) != 1:
                msg = ('I expected exactly one match for "CV Error" in '
                       'admixture\'s logfile, but I found {}. '
                       'Check: {}'.format(len(cv_error_lines), log_fn))
                raise Exception(msg)

            cv_error_match = re.search('(\d\.\d*$)', cv_error_lines[0])

        return float(cv_error_match.group(0))

    def _call_admixture(self, K, cores):
        command = '{} --cv {}.bed {} -j{}'
        command = command.format(self._EXECUTABLE, self.dataset.bedfile, K,
                                 cores)
        working_dir = getcwd()
        try:
            chdir(self.dataset.source.datasets_dir)
            with open(self.logfiles[K], 'w+') as logfile:
                subprocess.call(command.split(' '), stdout=logfile)
        finally:
            # Make sure you come back to the original working directory, even
            # in case of a KeyboardInterrupt or any Exception in general.
            chdir(working_dir)

    def _output_filepath(self, K, output_label):
        return self.dataset.bedfile + '.{}.{}'.format(K, output_label)
