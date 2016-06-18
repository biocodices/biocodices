import re
from glob import glob
from os.path import join, abspath, basename, expanduser
import pandas as pd
from termcolor import colored

from biocodices.components import Sample
from biocodices.variant_calling import VcfMunger
from biocodices.helpers import plural
from biocodices.plotters import AlignmentMetricsPlotter


class Cohort:
    def __init__(self, base_dir):
        """
        Expects the path to a directory that will have a 'data' subdirectory
        with fastq forward (R1) and reverse (R2) files from samples in it.
        """
        self.dir = abspath(expanduser(base_dir))
        self.id = basename(self.dir)
        self.data_dir = join(self.dir, 'data')
        self.results_dir = join(self.dir, 'results')
        self.samples = self._search_samples()
        self.sequencer_runs = list(set([sample.sequencer_run_id
                                        for sample in self.samples]))
        self.vcf_munger = VcfMunger()

        if len(self.samples) == 0:
            msg = 'I found no sample files (.fastq) in {}'
            raise EmptyCohortException(msg.format(self.data_dir))

        self.unfiltered_vcf = join(self.results_dir, 'raw_variants.vcf')
        self.filtered_vcf = join(self.results_dir, 'filtered.vcf')
        self.geno_filtered_vcf = join(self.results_dir,
                                      'filtered.geno-filtered.vcf')
        self.limited_vcf = join(self.results_dir,
                                'filtered.geno-filtered.lim.vcf')
        self.__vcf_stats = None

    def __repr__(self):
        tmpl = '{} with {} from {}'
        return tmpl.format(self.__class__.__name__,
                           plural('sample', len(self.samples)),
                           ', '.join(self.sequencer_runs))

    def vcf_stats(self, vcf_FORMAT_field=None):
        """Create a pandas DataFrame with the values from one FORMAT field of
        the passed VCF file: 'DP' for depth, 'GQ' for genotype quality.
        The df is memoized after the first call."""
        if self.__vcf_stats is None:
            _, samples_df = VcfMunger.vcf_to_frames(self.filtered_vcf)
            self.__vcf_stats = samples_df

        if vcf_FORMAT_field:
            col_index = pd.IndexSlice[:, vcf_FORMAT_field]
            return self.__vcf_stats.loc[:, col_index]

        return self.__vcf_stats

    def plot_alignment_metrics(self):
        frames = [sample.read_alignment_metrics() for sample in self.samples]
        df = pd.concat(frames, ignore_index=True).dropna(axis=1)
        # Remove the info of both pairs, leave by first and second of pair.
        df_pairs = df[df['CATEGORY'] != 'PAIR']
        plotter = AlignmentMetricsPlotter(df_pairs)
        plotter.plot_and_savefig(self.results_dir)

    def median_coverages(self):
        median_coverages = {}
        for sample in self.samples:
            median_coverages[sample.id] = sample.get_median_coverage()
        median_coverages = pd.Series(median_coverages)
        return median_coverages

    #  def plot_variant_calling_metrics(self):
        #  frames = [sample.read_variant_calling_metrics()
                  #  for sample in self.samples]
        #  df = pd.concat(frames, ignore_index=True).dropna(axis=1)
        #  plotter = VariantCallingMetricsPlotter(df)
        #  plotter.plot_and_savefig(self.dir)

    def _search_samples(self):
        glob_expr = join(self.data_dir, '*.{}'.format(Sample.reads_format))
        all_reads_filenames = sorted(glob(glob_expr))
        sample_ids = [re.search(r'(.*).R1', basename(reads_fn)).groups(1)[0]
                      for reads_fn in all_reads_filenames
                      if 'R1' in reads_fn]

        return [Sample(sample_id, self) for sample_id in sample_ids]

    def msg(self, msg):
        prefix = colored('[{}]'.format(self.id), 'magenta')
        return '{} {}'.format(prefix, msg)


class EmptyCohortException(Exception):
    pass
