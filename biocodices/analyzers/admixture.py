import re
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

from os import chdir, getcwd
from os.path import expanduser, isfile
from biocodices.plotters.admixture_plotter import AdmixturePlotter
from biocodices.helpers.config import Config


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
        self.plotter = AdmixturePlotter(self, self.dataset.dir)

    def __repr__(self):
        return "<Admixture for {}, Ks={}>".format(self.dataset.label,
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

            #  regions = self.result[K].index.get_level_values('region').unique()
            #  if infer_components and len(regions) >= 3:
                #  self._assign_regions_to_clusters(self.result[K])
            #  self._infer_clusters_from_reference_population(self.result[K])
            # self.result[K] = self._reorder_clusters(self.result[K])

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

    def _read_ancestry_file(self, K):
        result = pd.read_table(self.Qfiles[K], sep='\s+', header=None)
        result.columns = range(1, K+1)
        result.index = self.dataset.samplegroup.samples.index
        result = result.join(self.dataset.samplegroup.samples)
        multi_index = ['sample']
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
        command = '{} --cv {} {} -j{}'
        command = command.format(self._EXECUTABLE, self.dataset.bed, K,
                                 cores)
        working_dir = getcwd()
        try:
            chdir(self.dataset.dir)
            with open(self.logfiles[K], 'w+') as logfile:
                subprocess.call(command.split(' '), stdout=logfile)
        finally:
            # Make sure you come back to the original working directory, even
            # in case of a KeyboardInterrupt or any Exception in general.
            chdir(working_dir)

    def _output_filepath(self, K, output_label):
        return self.dataset.path_label + '.{}.{}'.format(K, output_label)
